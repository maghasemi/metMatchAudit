# Checks each row against the RefMet entry its matching step assigned it
# to (notes/decisions.md, 2026-10-07, v0.5 plan). MetLinkR's standardized
# name is a RefMet name, so every group stands for one RefMet entry. The
# check asks whether the row's own evidence agrees with that entry. It
# reports, and does not change any pair result.

#' Load the RefMet table
#'
#' Reads the text download of the full RefMet table
#' (`https://www.metabolomicsworkbench.org/rest/refmet/all/xx/all/txt`):
#' one `field<TAB>value` line per field, a blank line between entries,
#' each entry starting with `name`.
#'
#' @param path Path to the downloaded file.
#' @param expected_md5 If given, the run stops when the file's md5
#'   differs, so a new RefMet release is never picked up by accident.
#'
#' @return A data frame with one row per entry and the columns `name`,
#'   `refmet_id`, `formula`, `inchi_key`, `pubchem_cid`, `super_class`,
#'   `main_class`, `sub_class`, with class `refmet_db` and the file's md5
#'   in attribute `md5`. Invalid InChIKeys become `NA`.
#'
#' @export
load_refmet <- function(path, expected_md5 = NULL) {
  if (!file.exists(path)) stop("Cannot find ", path, ".", call. = FALSE)
  md5 <- unname(tools::md5sum(path))
  if (!is.null(expected_md5) && md5 != expected_md5) {
    stop("md5 of ", path, " is ", md5, ", expected ", expected_md5, ".", call. = FALSE)
  }
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- x[nzchar(x)]
  tab <- regexpr("\t", x, fixed = TRUE)
  field <- ifelse(tab > 0, substr(x, 1, tab - 1), x)
  value <- ifelse(tab > 0, substr(x, tab + 1, nchar(x)), NA_character_)
  value <- trimws(value)
  value[!is.na(value) & value == ""] <- NA_character_
  rec <- cumsum(field == "name")
  if (rec[1] == 0) stop(path, " does not start with a `name` line.", call. = FALSE)

  cols <- c("name", "refmet_id", "formula", "inchi_key", "pubchem_cid",
            "super_class", "main_class", "sub_class")
  n <- max(rec)
  out <- data.frame(matrix(NA_character_, n, length(cols), dimnames = list(NULL, cols)),
                    stringsAsFactors = FALSE)
  for (col in cols) {
    sel <- field == col
    out[[col]][rec[sel]] <- value[sel]
  }
  out$inchi_key[!is_valid_inchikey(out$inchi_key)] <- NA_character_
  if (anyDuplicated(out$name)) {
    stop("Some RefMet names are used by more than one entry.", call. = FALSE)
  }
  class(out) <- c("refmet_db", "data.frame")
  attr(out, "md5") <- md5
  out
}

# Parent compounds (salts and counterions removed) with their InChIKey and
# formula. Uses the PubChem client and cache from identifiers.R.
.props_for_cids <- function(cids, cache, fetch, batch = 100) {
  cids <- unique(cids[!is.na(cids)])
  out <- data.frame(cid = cids, key = NA_character_, formula = NA_character_,
                    stringsAsFactors = FALSE)
  if (length(cids) == 0) return(out)
  for (start in seq(1, length(cids), by = batch)) {
    chunk <- cids[start:min(start + batch - 1, length(cids))]
    res <- .cached(cache, fetch, "compound/cid/property/InChIKey,MolecularFormula/JSON",
                   body = list(cid = paste(chunk, collapse = ",")))
    props <- res$PropertyTable$Properties
    if (!is.null(props) && nrow(props) > 0) {
      i <- match(as.character(props$CID), out$cid)
      out$key[i[!is.na(i)]] <- props$InChIKey[!is.na(i)]
      out$formula[i[!is.na(i)]] <- props$MolecularFormula[!is.na(i)]
    }
  }
  out
}

# One parent for a key or a CID, or NA when PubChem has none or gives
# parents that disagree on block1.
.one_parent <- function(parent_cids, cache, fetch) {
  none <- list(key = NA_character_, formula = NA_character_)
  if (length(parent_cids) == 0) return(none)
  p <- .props_for_cids(parent_cids, cache, fetch)
  p <- p[!is.na(p$key), , drop = FALSE]
  if (nrow(p) == 0 || length(unique(substr(p$key, 1, 14))) != 1) return(none)
  list(key = p$key[1], formula = if (length(unique(p$formula)) == 1) p$formula[1] else NA_character_)
}

.parent_of_key <- function(key, cache, fetch) {
  res <- .cached(cache, fetch, paste0("compound/inchikey/", key, "/cids/JSON?cids_type=parent"))
  cids <- res$IdentifierList$CID
  .one_parent(if (is.null(cids)) character(0) else as.character(cids), cache, fetch)
}

.parent_of_cid <- function(cid, cache, fetch) {
  res <- .cached(cache, fetch, paste0("compound/cid/", cid, "/cids/JSON?cids_type=parent"))
  cids <- res$IdentifierList$CID
  .one_parent(if (is.null(cids)) character(0) else as.character(cids), cache, fetch)
}

#' Check each row against its RefMet entry
#'
#' The group name of a row is taken as the RefMet entry the matching step
#' assigned it to. The first case that applies decides how the row is
#' compared with that entry:
#'
#' 1. The entry is a lipid name that does not fix one structure: rgoslin
#'    reads it at species, molecular species or sn-position level, such as
#'    `"CAR 5:0"` or `"MG 18:2"`. Its InChIKey stands for one member and
#'    is not used, for the same reason [resolve_identifiers()] skips such
#'    names. A row with a lipid profile is compared by Species. A row with
#'    only a key is compared by molecular formula: the formula of the
#'    row's parent compound against the entry's formula.
#' 2. The entry is a lipid at structure level and the row has a lipid
#'    profile: compared by Species.
#' 3. The entry has an InChIKey and the row has a key: compared by
#'    block1. When the raw keys differ, both are first taken to their
#'    PubChem parent compounds, which removes salts and counterions.
#'
#' Species that differ only in hexose headgroup naming (HexCer against
#' GlcCer) are not a conflict, as in rule R04h. A difference becomes
#' `"identifier_conflict"` when the row's name and the entry's name are the
#' same after [normalize_name()], and `"conflict"` otherwise.
#'
#' Rows with group names, alternative names or internal standards are
#' skipped.
#'
#' @param evidence Rows after [row_evidence()] and [lipid_evidence()],
#'   with a `group_id` column.
#' @param refmet A table from [load_refmet()].
#' @param cache_file,fetch As in [resolve_identifiers()].
#'
#' @return `evidence` with columns added:
#'   * `rm_check`: `"agrees"`, `"conflict"`, `"identifier_conflict"`,
#'     `"unchecked"` (a needed PubChem answer was missing or ambiguous),
#'     `"not_compared"`, `"skipped"`, `"entry_missing"` or `"no_group"`
#'   * `rm_basis`: `"species"`, `"formula"`, `"block1"` or
#'     `"parent_block1"`
#'   * `rm_detail`: the values compared
#'
#' @export
refmet_check <- function(evidence, refmet, cache_file = NULL, fetch = NULL) {
  if (!inherits(refmet, "refmet_db")) stop("`refmet` must come from load_refmet().", call. = FALSE)
  needed <- c("row_id", "original_name", "group_id", "inchikey", "ik_valid", "ik_block1",
              "lipid_species", "name_group", "name_alternatives", "name_internal_standard")
  missing <- setdiff(needed, names(evidence))
  if (length(missing) > 0) {
    stop("`evidence` is missing ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  n <- nrow(evidence)
  check  <- rep(NA_character_, n)
  basis  <- rep(NA_character_, n)
  detail <- rep(NA_character_, n)

  g <- as.character(evidence$group_id)
  check[is.na(g) | g %in% c("", "-")] <- "no_group"
  e <- match(g, refmet$name)
  check[is.na(check) & is.na(e)] <- "entry_missing"
  skip <- evidence$name_group | evidence$name_alternatives | evidence$name_internal_standard
  check[is.na(check) & skip] <- "skipped"

  todo <- which(is.na(check))
  entry_names <- unique(refmet$name[e[todo]])
  lp <- parse_lipid_names(entry_names)
  lp_level <- ifelse(lp$parsed, lp$level, NA_character_)
  lp_species <- ifelse(lp$parsed, lp$species, NA_character_)
  ent_level <- rep(NA_character_, n); ent_species <- rep(NA_character_, n)
  ent_level[todo]   <- lp_level[match(refmet$name[e[todo]], entry_names)]
  ent_species[todo] <- lp_species[match(refmet$name[e[todo]], entry_names)]

  row_sp  <- evidence$lipid_species
  has_key <- !is.na(evidence$ik_valid) & evidence$ik_valid
  ent_key <- rep(NA_character_, n); ent_key[todo] <- refmet$inchi_key[e[todo]]

  case <- rep(NA_character_, n)
  # A lipid name below structure level does not fix one structure, so the
  # entry's key is only one member (same reasoning as in identifiers.R).
  sp_entry <- !is.na(ent_level) & ent_level %in% c("species", "molecular", "sn_position")
  fine_entry <- !is.na(ent_level) & ent_level == "structure"
  case[todo] <- ifelse(sp_entry[todo] & !is.na(row_sp[todo]), "species",
                ifelse(sp_entry[todo] & has_key[todo], "formula",
                ifelse(sp_entry[todo], "none",
                ifelse(fine_entry[todo] & !is.na(row_sp[todo]), "species",
                ifelse(!is.na(ent_key[todo]) & has_key[todo], "key", "none")))))
  check[!is.na(case) & case == "none"] <- "not_compared"

  differs <- rep(NA, n)

  sp <- which(case == "species")
  if (length(sp) > 0) {
    same <- row_sp[sp] == ent_species[sp] |
      .generic_headgroup(row_sp[sp]) == .generic_headgroup(ent_species[sp])
    differs[sp] <- !same
    basis[sp] <- "species"
    detail[sp] <- sprintf("row Species %s, entry Species %s", row_sp[sp], ent_species[sp])
  }

  # Rows whose raw key already shares block1 with the entry need no lookup.
  ky <- which(case == "key")
  raw_same <- ky[substr(ent_key[ky], 1, 14) == evidence$ik_block1[ky]]
  differs[raw_same] <- FALSE
  basis[raw_same] <- "block1"
  detail[raw_same] <- sprintf("row key %s, entry key %s", evidence$inchikey[raw_same], ent_key[raw_same])

  need_parent <- c(setdiff(ky, raw_same), which(case == "formula"))
  if (length(need_parent) > 0) {
    if (is.null(fetch)) {
      if (!requireNamespace("httr2", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
        message("httr2 and jsonlite are needed for parent compounds; ",
                length(need_parent), " rows stay unchecked.")
        fetch <- function(path, body = NULL) list(error = "no client")
      } else {
        fetch <- .pubchem_fetch
      }
    }
    cache <- .new_cache(cache_file)
    for (i in need_parent) {
      rp <- .parent_of_key(evidence$inchikey[i], cache, fetch)
      if (case[i] == "formula") {
        ef <- refmet$formula[e[i]]
        basis[i] <- "formula"
        if (is.na(rp$formula) || is.na(ef)) { check[i] <- "unchecked"; next }
        differs[i] <- rp$formula != ef
        detail[i] <- sprintf("row parent formula %s, entry formula %s", rp$formula, ef)
      } else {
        cid <- refmet$pubchem_cid[e[i]]
        ep <- if (!is.na(cid)) .parent_of_cid(cid, cache, fetch) else .parent_of_key(ent_key[i], cache, fetch)
        basis[i] <- "parent_block1"
        if (is.na(rp$key) || is.na(ep$key)) { check[i] <- "unchecked"; next }
        differs[i] <- substr(rp$key, 1, 14) != substr(ep$key, 1, 14)
        detail[i] <- sprintf("row key %s (parent %s), entry key %s (parent %s)",
                             evidence$inchikey[i], rp$key, ent_key[i], ep$key)
      }
    }
  }

  done <- which(!is.na(differs) & is.na(check))
  same_name <- normalize_name(evidence$original_name[done]) == normalize_name(refmet$name[e[done]])
  check[done] <- ifelse(!differs[done], "agrees",
                 ifelse(!is.na(same_name) & same_name, "identifier_conflict", "conflict"))

  evidence$rm_entry  <- ifelse(is.na(e), NA_character_, refmet$name[e])
  evidence$rm_check  <- check
  evidence$rm_basis  <- basis
  evidence$rm_detail <- detail
  evidence
}
