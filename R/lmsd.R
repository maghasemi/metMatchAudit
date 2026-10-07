# LIPID MAPS Structure Database (LMSD), a second lipid source looked up
# by InChIKey the same way as SwissLipids. Tested in
# scripts/82_lmsd_feasibility.R (notes/decisions.md, 2026-10-07): where
# SwissLipids and LMSD both gave a Species, they agreed in 149 of 149
# development rows and 66 of 66 COMETS rows.

#' Load the LIPID MAPS Structure Database
#'
#' Reads the standard LMSD structure-data file (`structures.sdf` inside
#' the zip from <https://www.lipidmaps.org/databases/lmsd/download>) and
#' keeps, for every record with a valid InChIKey and an abbreviation, the
#' abbreviation under its full key, its first 25 characters and its first
#' 14. Reading the file takes about a minute. The object it returns can
#' be kept with [saveRDS()].
#'
#' LMSD is released under CC BY 4.0. Cite LIPID MAPS when you publish
#' results that use it.
#'
#' @param path Path to the SDF file.
#' @param expected_md5 Optional md5 checksum. If given and the file does
#'   not match, loading stops.
#'
#' @return An object of class `lmsd_db`.
#'
#' @export
load_lmsd <- function(path, expected_md5 = NULL) {
  if (!file.exists(path)) stop("Cannot find ", path, ".", call. = FALSE)
  md5 <- unname(tools::md5sum(path))
  if (!is.null(expected_md5) && !identical(md5, expected_md5)) {
    stop("md5 of ", path, " is ", md5, ", expected ", expected_md5, ".", call. = FALSE)
  }
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  # Each field is a "> <NAME>" line followed by its value; "$$$$" ends a record.
  rec <- cumsum(c(0L, utils::head(x == "$$$$", -1))) + 1L
  field_value <- function(f) {
    out <- rep(NA_character_, max(rec))
    i <- which(x == paste0("> <", f, ">"))
    out[rec[i]] <- trimws(x[i + 1])
    out
  }
  key  <- field_value("INCHI_KEY")
  abbr <- field_value("ABBREVIATION")
  n_records <- sum(x == "$$$$")
  rm(x)
  ok <- is_valid_inchikey(key) & !is.na(abbr) & nzchar(abbr)
  key <- key[ok]; abbr <- abbr[ok]

  db <- list(
    md5       = md5,
    n_records = n_records,
    n_used    = length(key),
    by_key    = .as_env_list(key, abbr),
    by_key25  = .as_env_list(substr(key, 1, 25), abbr),
    by_key14  = .as_env_list(substr(key, 1, 14), abbr)
  )
  class(db) <- "lmsd_db"
  db
}

#' @export
print.lmsd_db <- function(x, ...) {
  cat("<lmsd_db> ", x$n_records, " records, ", x$n_used,
      " with a key and an abbreviation, md5 ", x$md5, "\n", sep = "")
  invisible(x)
}

#' Look up InChIKeys in LMSD
#'
#' Tries the full key, then its first 25 characters, then its first 14,
#' and stops at the first tier that finds anything. The abbreviations of
#' all hits are parsed with [parse_lipid_names()]. A key gets a Species
#' only when every hit parses and all give the same Species, and a
#' Molecular species only when they also agree on that.
#'
#' @param inchikeys A character vector of InChIKeys.
#' @param db An object from [load_lmsd()].
#'
#' @return A data frame with one row per key: `lm_tier` (1, 2 or 3, or
#'   `NA` without a hit), `lm_n_hits`, `lm_species`, `lm_molecular`.
#'
#' @export
lmsd_lookup <- function(inchikeys, db) {
  if (!inherits(db, "lmsd_db")) stop("`db` must come from load_lmsd().", call. = FALSE)
  inchikeys <- as.character(inchikeys)
  hits <- lapply(inchikeys, function(k) {
    if (is.na(k) || !is_valid_inchikey(k)) return(list(tier = NA_integer_, abbr = character(0)))
    for (t in 1:3) {
      env <- switch(t, db$by_key, db$by_key25, db$by_key14)
      a <- .get_many(switch(t, k, substr(k, 1, 25), substr(k, 1, 14)), env)
      a <- unique(a[!is.na(a)])
      if (length(a) > 0) return(list(tier = t, abbr = a))
    }
    list(tier = NA_integer_, abbr = character(0))
  })
  all_abbr <- unique(unlist(lapply(hits, `[[`, "abbr")))
  lp <- parse_lipid_names(all_abbr)
  sp  <- stats::setNames(ifelse(lp$parsed, lp$species, NA_character_), all_abbr)
  mol <- stats::setNames(ifelse(lp$parsed & lp$level != "species", lp$molecular, NA_character_), all_abbr)
  one <- function(v) { v <- unique(v); if (length(v) == 1 && !is.na(v)) v else NA_character_ }
  species <- vapply(hits, function(h) if (length(h$abbr) == 0) NA_character_ else one(sp[h$abbr]), character(1))
  data.frame(
    lm_tier      = vapply(hits, `[[`, integer(1), "tier"),
    lm_n_hits    = vapply(hits, function(h) length(h$abbr), integer(1)),
    lm_species   = species,
    lm_molecular = ifelse(is.na(species), NA_character_,
                          vapply(hits, function(h) if (length(h$abbr) == 0) NA_character_ else one(mol[h$abbr]),
                                 character(1))),
    stringsAsFactors = FALSE
  )
}
