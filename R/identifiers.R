# Identifier resolution: PubChem CIDs and HMDB IDs to InChIKeys.
#
# Tested in scripts/07_v0.4/02_id_resolution_feasibility.R and scripts/07_v0.4/03_id_resolution_retest.R (notes/decisions.md, 2026-10-06):
# on 783 rows not used in development, CID-derived keys agreed with the
# stored key in 779 of 779 cases and HMDB-derived keys in 288 of 297.
# A derived key describes the database record for the identifier. If a
# vendor put the wrong identifier in a row, the error is passed on, which
# is why every key keeps its source.

.id_tokens <- function(x) {
  lapply(as.character(x), function(cell) {
    if (is.na(cell) || !nzchar(trimws(cell))) return(character(0))
    parts <- trimws(unlist(strsplit(cell, "[,;| ]+")))
    parts <- sub("^(pubchem|hmdb|cid)[:_]", "", parts, ignore.case = TRUE)
    parts[nzchar(parts)]
  })
}

#' Normalise PubChem CID cells
#'
#' Removes prefixes such as `"pubchem:"`. A cell that lists several IDs
#' gives `NA`, because it does not say which compound the row is.
#'
#' @param x A character vector.
#' @return A character vector of single CIDs or `NA`.
#' @examples
#' normalize_pubchem_id(c("pubchem:96", "12594", "pubchem:439685,7405", NA))
#' @export
normalize_pubchem_id <- function(x) {
  vapply(.id_tokens(x), function(t) {
    t <- unique(t)
    if (length(t) == 1 && grepl("^[0-9]+$", t)) t else NA_character_
  }, character(1))
}

#' Normalise HMDB ID cells
#'
#' Removes prefixes such as `"hmdb:"` and writes the ID in the current
#' seven-digit form: `"HMDB00641"` and `"HMDB641"` both become
#' `"HMDB0000641"`. A cell listing several IDs gives `NA`.
#'
#' @param x A character vector.
#' @return A character vector of single HMDB IDs or `NA`.
#' @examples
#' normalize_hmdb_id(c("HMDB00641", "hmdb:HMDB14493", "hmdb:HMDB0240679,HMDB0240680"))
#' @export
normalize_hmdb_id <- function(x) {
  vapply(.id_tokens(x), function(t) {
    t <- unique(toupper(t))
    if (length(t) != 1 || !grepl("^HMDB[0-9]+$", t)) return(NA_character_)
    digits <- sub("^HMDB", "", t)
    if (nchar(digits) > 7) return(NA_character_)
    paste0("HMDB", formatC(as.integer(digits), width = 7, flag = "0"))
  }, character(1))
}

# Default way of asking PubChem. Returns parsed JSON, or a list with an
# `error` element.
.pubchem_fetch <- function(path, body = NULL) {
  req <- httr2::request(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/", path))
  req <- httr2::req_user_agent(req, "metMatchAudit")
  req <- httr2::req_timeout(req, 60)
  req <- httr2::req_error(req, is_error = function(r) FALSE)
  if (!is.null(body)) req <- do.call(httr2::req_body_form, c(list(req), body))
  resp <- tryCatch(httr2::req_perform(req), error = function(e) e)
  Sys.sleep(0.25)  # PubChem asks for at most 5 requests per second
  if (inherits(resp, "error")) return(list(error = conditionMessage(resp)))
  tryCatch(jsonlite::fromJSON(httr2::resp_body_string(resp)),
           error = function(e) list(error = "response was not JSON"))
}

# A request has failed when it gives an error, or a PubChem Fault other
# than "not found", such as PUGREST.ServerBusy. Not found is a real answer.
.failed <- function(res) {
  if (!is.null(res$error)) return(TRUE)
  if (is.null(res$Fault)) return(FALSE)
  !identical(res$Fault$Code, "PUGREST.NotFound")
}

# A cache is an environment, written to `file` after each new answer so
# an interrupted run loses nothing. `failed` counts the requests of this
# run that failed. Their answers are not kept, so the next run asks again.
.new_cache <- function(file) {
  env <- if (!is.null(file) && file.exists(file)) readRDS(file) else new.env()
  failed <- new.env()
  failed$n <- 0
  list(env = env, file = file, failed = failed)
}

.cached <- function(cache, fetch, path, body = NULL) {
  key <- paste(path, paste(unlist(body), collapse = ","))
  hit <- cache$env[[key]]
  # A failed answer stored by a version before 0.9.1 is asked again.
  if (!is.null(hit) && !.failed(hit)) return(hit)
  res <- fetch(path, body)
  if (.failed(res)) {
    cache$failed$n <- cache$failed$n + 1
  } else {
    cache$env[[key]] <- res
    if (!is.null(cache$file)) {
      dir.create(dirname(cache$file), showWarnings = FALSE, recursive = TRUE)
      saveRDS(cache$env, cache$file)
    }
  }
  res
}

.parent_cid <- function(cid, cache, fetch) {
  res <- .cached(cache, fetch, paste0("compound/cid/", cid, "/cids/JSON?cids_type=parent"))
  # Without an answer the parent is unknown. Using the CID itself would
  # give the key of what may be a salt, and it would look like an answer.
  if (.failed(res)) return(NA_character_)
  p <- res$IdentifierList$CID
  if (is.null(p) || length(p) == 0) cid else as.character(p[1])
}

.keys_for_cids <- function(cids, cache, fetch, batch = 100) {
  cids <- unique(cids[!is.na(cids)])
  out <- stats::setNames(rep(NA_character_, length(cids)), cids)
  if (length(cids) == 0) return(out)
  for (start in seq(1, length(cids), by = batch)) {
    chunk <- cids[start:min(start + batch - 1, length(cids))]
    res <- .cached(cache, fetch, "compound/cid/property/InChIKey/JSON",
                   body = list(cid = paste(chunk, collapse = ",")))
    props <- res$PropertyTable$Properties
    if (!is.null(props) && nrow(props) > 0) out[as.character(props$CID)] <- props$InChIKey
  }
  out
}

.hmdb_cids <- function(hmdb, cache, fetch) {
  res <- .cached(cache, fetch, paste0("compound/xref/RegistryID/", hmdb, "/cids/JSON"))
  cids <- res$IdentifierList$CID
  if (is.null(cids)) character(0) else as.character(cids)
}

#' Fill in missing InChIKeys from PubChem CIDs or HMDB IDs
#'
#' For each row without a valid InChIKey, tries its PubChem CID first and
#' its HMDB ID second:
#'
#' * CID: the InChIKey of PubChem's parent compound for that CID (salts
#'   and counterions removed, neutralised).
#' * HMDB: the CIDs PubChem lists for the HMDB registry ID, each taken to
#'   its parent compound. Used only if they all give the same key.
#'
#' Rows with a group name or alternative names are not resolved, and
#' neither are cells that list several IDs. Lipid names that rgoslin
#' reads at Species, molecular species or sn-position level (such as
#' "TG 42:0" or "PC 34:1") are not resolved either: no single structure
#' matches them, so a key from their identifier would claim more than the
#' row reports. Internal standards are skipped too, since their
#' identifier is usually that of the unlabelled compound. Rows that
#' already have a valid key keep it.
#'
#' Needs httr2, jsonlite and internet access, unless every answer is
#' already in the cache.
#'
#' A request that fails (a network error, or a PubChem error such as a
#' busy server) is not cached, and the rows that needed it get no key in
#' this run. A message gives the number of failed requests, and running
#' again asks them again. "Not found" is a real answer and is cached.
#'
#' @param rows A data frame with `original_name` and `inchikey`, and
#'   optionally `pubchem` and `hmdb` columns.
#' @param cache_file Where to keep PubChem answers between runs. `NULL`
#'   keeps them only for this call.
#' @param fetch The function that asks PubChem. Only changed in tests.
#'
#' @return `rows` with `inchikey` filled in where possible, the original
#'   value kept in `inchikey_file`, and `inchikey_source`: `"file"`,
#'   `"pubchem_cid"`, `"hmdb_via_pubchem"`, or `NA`. The number of failed
#'   PubChem requests is in `attr(, "pubchem_failed")`; it is absent when
#'   no row needed PubChem.
#'
#' @export
resolve_identifiers <- function(rows, cache_file = NULL, fetch = NULL) {
  for (col in c("original_name", "inchikey")) {
    if (!col %in% names(rows)) stop("`rows` needs a column called `", col, "`.", call. = FALSE)
  }
  if ("inchikey_source" %in% names(rows)) {
    stop("`rows` already has `inchikey_source`; identifiers were resolved before.", call. = FALSE)
  }
  n <- nrow(rows)
  rows$inchikey_file <- rows$inchikey
  has_key <- is_valid_inchikey(rows$inchikey)
  rows$inchikey_source <- ifelse(has_key, "file", NA_character_)

  cid  <- if ("pubchem" %in% names(rows)) normalize_pubchem_id(rows$pubchem) else rep(NA_character_, n)
  hmdb <- if ("hmdb" %in% names(rows)) normalize_hmdb_id(rows$hmdb) else rep(NA_character_, n)
  # A lipid name that rgoslin reads as a Species, molecular species or
  # sn-position annotation does not describe one structure. An identifier
  # on such a row points to one representative record, and its key would
  # claim a structure the row never reported. Found on the development
  # data: Broad rows such as "TG 42:0" or "SM 34:2;O2" got keys that
  # conflicted with their own InChI and with their lipid profile.
  lipid_level <- parse_lipid_names(rows$original_name)$level
  not_one_structure <- !is.na(lipid_level) & lipid_level %in% c("species", "molecular", "sn_position")
  # Internal standards carry the identifier of the unlabelled compound,
  # so a key from it would describe the wrong molecule (found on Broad
  # "Inosine-15N4 [iSTD]", which got the key of plain inosine).
  names_ok <- !is_group_name(rows$original_name) & !has_alternative_names(rows$original_name) &
    !not_one_structure & !is_internal_standard(rows$original_name)
  todo_cid  <- which(!has_key & names_ok & !is.na(cid))
  todo_hmdb <- which(!has_key & names_ok & is.na(cid) & !is.na(hmdb))
  if (length(todo_cid) + length(todo_hmdb) == 0) return(rows)

  if (is.null(fetch)) {
    if (!requireNamespace("httr2", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
      message("httr2 and jsonlite are needed to resolve identifiers; ",
              length(todo_cid) + length(todo_hmdb), " rows stay without a key.")
      return(rows)
    }
    fetch <- .pubchem_fetch
  }
  cache <- .new_cache(cache_file)

  if (length(todo_cid) > 0) {
    parents <- vapply(cid[todo_cid], .parent_cid, character(1), cache = cache, fetch = fetch)
    keys <- .keys_for_cids(parents, cache, fetch)
    found <- unname(keys[parents])
    ok <- !is.na(found)
    rows$inchikey[todo_cid[ok]] <- found[ok]
    rows$inchikey_source[todo_cid[ok]] <- "pubchem_cid"
  }

  # Rows whose CID did not resolve also get a chance through HMDB.
  todo_hmdb <- union(todo_hmdb, todo_cid[is.na(rows$inchikey_source[todo_cid]) & !is.na(hmdb[todo_cid])])
  for (i in todo_hmdb) {
    failed_before <- cache$failed$n
    cids <- .hmdb_cids(hmdb[i], cache, fetch)
    if (length(cids) == 0) next
    parents <- vapply(cids, .parent_cid, character(1), cache = cache, fetch = fetch)
    k <- unique(stats::na.omit(unname(.keys_for_cids(parents, cache, fetch)[parents])))
    # All CIDs must agree, so one that could not be asked leaves the row open.
    if (cache$failed$n > failed_before) next
    if (length(k) == 1) {
      rows$inchikey[i] <- k
      rows$inchikey_source[i] <- "hmdb_via_pubchem"
    }
  }

  n_failed <- cache$failed$n
  if (n_failed > 0) {
    message(n_failed, " PubChem request", if (n_failed > 1) "s" else "", " failed. ",
            "Rows that needed ", if (n_failed > 1) "them" else "it", " have no key from PubChem in this run. ",
            "Run again to ask again", if (is.null(cache_file)) "; give `cache_file` to keep the answers already received" else "",
            ".")
  }
  attr(rows, "pubchem_failed") <- n_failed
  rows
}
