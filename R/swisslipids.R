# SwissLipids, read from its downloadable lipids.tsv file.
#
# The file is large (about 690 MB, 780k rows). load_swisslipids() reads
# it once and builds lookup tables, so later calls are fast. The object
# it returns can be kept with saveRDS() between sessions.

# Columns we use, by position in lipids.tsv, and the header text
# expected at each position. Names are compared after trimming, since
# two headers in the file end with a space.
.sl_columns <- c(
  "1"  = "Lipid ID",
  "2"  = "Level",
  "4"  = "Abbreviation*",
  "7"  = "Parent",
  "11" = "InChI key (pH7.3)"
)

#' Load the SwissLipids table
#'
#' Reads `lipids.tsv` as downloaded from
#' <https://www.swisslipids.org/#/downloads> and prepares it for lookups
#' by InChIKey and for walking up the lipid hierarchy.
#'
#' SwissLipids is released under CC BY 4.0. Cite it when you publish
#' results that use it.
#'
#' @param path Path to `lipids.tsv`.
#' @param expected_md5 Optional md5 checksum. If given and the file does
#'   not match, loading stops. Use it to make sure an analysis runs on
#'   the same release it was built on.
#'
#' @return An object of class `swisslipids_db`. Its `md5` element records
#'   the checksum of the file that was read.
#'
#' @export
load_swisslipids <- function(path, expected_md5 = NULL) {
  if (!file.exists(path)) {
    stop("Cannot find ", path, ".", call. = FALSE)
  }
  md5 <- unname(tools::md5sum(path))
  if (!is.null(expected_md5) && !identical(md5, expected_md5)) {
    stop("md5 of ", path, " is ", md5, ", expected ", expected_md5, ".",
         call. = FALSE)
  }

  header <- strsplit(readLines(path, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  pos <- as.integer(names(.sl_columns))
  found <- trimws(header[pos])
  if (any(is.na(found)) || !identical(found, unname(.sl_columns))) {
    stop("The header of ", path, " does not match the SwissLipids layout ",
         "this package expects. Expected ",
         paste(sprintf("'%s' at column %d", .sl_columns, pos), collapse = ", "),
         ".", call. = FALSE)
  }

  # quote = "" because names and SMILES can contain quote characters.
  tab <- data.table::fread(path, sep = "\t", quote = "", select = pos,
                           colClasses = "character", na.strings = "",
                           encoding = "UTF-8", showProgress = FALSE,
                           data.table = FALSE)
  names(tab) <- c("slm_id", "level", "abbrev", "parent", "inchikey")

  tab$inchikey <- sub("^InChIKey=", "", tab$inchikey)
  tab$inchikey[!is_valid_inchikey(tab$inchikey)] <- NA_character_

  with_key <- !is.na(tab$inchikey)
  ids  <- tab$slm_id[with_key]
  keys <- tab$inchikey[with_key]

  db <- list(
    md5        = md5,
    n_rows     = nrow(tab),
    parent     = .as_env(tab$slm_id, tab$parent),
    level      = .as_env(tab$slm_id, tab$level),
    abbrev     = .as_env(tab$slm_id, tab$abbrev),
    by_key     = .as_env_list(keys, ids),
    by_key25   = .as_env_list(substr(keys, 1, 25), ids),
    by_key14   = .as_env_list(substr(keys, 1, 14), ids)
  )
  class(db) <- "swisslipids_db"
  db
}

#' @export
print.swisslipids_db <- function(x, ...) {
  cat("<swisslipids_db>", x$n_rows, "entries, md5", x$md5, "\n")
  invisible(x)
}

# One value per name.
.as_env <- function(keys, values) {
  keep <- !is.na(keys)
  list2env(as.list(stats::setNames(values[keep], keys[keep])), hash = TRUE)
}

# Several values per name, for keys that more than one entry shares.
.as_env_list <- function(keys, values) {
  list2env(split(values, keys), hash = TRUE)
}

.get_one <- function(ids, env) {
  if (length(ids) == 0) return(character(0))
  unlist(mget(ids, envir = env, ifnotfound = list(NA_character_)),
         use.names = FALSE)
}

.get_many <- function(key, env) {
  if (is.na(key)) return(character(0))
  out <- mget(key, envir = env, ifnotfound = list(character(0)))[[1]]
  out
}

.split_ids <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  ids <- trimws(unlist(strsplit(x, "[|,]")))
  ids[grepl("^SLM:", ids)]
}

# All ancestors of one entry, the entry itself first. At most 15 steps,
# and no entry is visited twice, so a loop in the parent links cannot
# hang.
.sl_lineage <- function(id, db) {
  seen <- id
  frontier <- id
  for (step in 1:15) {
    nxt <- setdiff(.split_ids(.get_one(frontier, db$parent)), seen)
    if (length(nxt) == 0) break
    seen <- c(seen, nxt)
    frontier <- nxt
  }
  seen
}

# Abbreviations of all entries at one level in the lineage of the
# given ids, sorted and without duplicates.
.sl_at_level <- function(ids, level, db) {
  out <- unlist(lapply(ids, function(id) {
    lin <- .sl_lineage(id, db)
    .get_one(lin[.get_one(lin, db$level) %in% level], db$abbrev)
  }))
  sort(unique(out[!is.na(out)]))
}

#' Place InChIKeys in the SwissLipids hierarchy
#'
#' For each key, finds matching SwissLipids entries and collects the
#' Species and Molecular subspecies above them.
#'
#' Matching runs in three tiers and stops at the first that finds
#' anything: the full key; the first 25 characters, which ignores the
#' protonation letter (SwissLipids stores structures as charged at
#' pH 7.3); the first 14 characters, which is the connectivity block.
#'
#' If none of the matched entries leads to a Species, all entries with
#' the same first 14 characters are tried instead. Same connectivity
#' means the same chain composition, so those entries should share one
#' Species. The fallback is used only if they name exactly one.
#'
#' @param inchikeys A character vector of InChIKeys. Invalid keys give an
#'   empty result row.
#' @param db An object from [load_swisslipids()].
#'
#' @return A data frame with one row per key: `sl_tier` (1, 2, 3, or
#'   `NA` if nothing matched), `sl_n_hits`, `sl_species` and
#'   `sl_molecular` (SwissLipids abbreviations, several joined by `"; "`),
#'   `sl_n_species`, and `sl_fallback` (`TRUE` if the Species came from
#'   the connectivity fallback).
#'
#' @export
swisslipids_lookup <- function(inchikeys, db) {
  if (!inherits(db, "swisslipids_db")) {
    stop("`db` must come from load_swisslipids().", call. = FALSE)
  }
  inchikeys <- as.character(inchikeys)
  rows <- lapply(inchikeys, function(k) .sl_one_key(k, db))
  do.call(rbind, rows)
}

.sl_one_key <- function(key, db) {
  empty <- data.frame(sl_tier = NA_integer_, sl_n_hits = 0L,
                      sl_species = NA_character_, sl_molecular = NA_character_,
                      sl_n_species = 0L, sl_fallback = FALSE,
                      stringsAsFactors = FALSE)
  if (!is_valid_inchikey(key)) return(empty)

  tiers <- list(.get_many(key, db$by_key),
                .get_many(substr(key, 1, 25), db$by_key25),
                .get_many(substr(key, 1, 14), db$by_key14))
  tier <- which(lengths(tiers) > 0)[1]
  if (is.na(tier)) return(empty)
  hits <- tiers[[tier]]

  species   <- .sl_at_level(hits, "Species", db)
  molecular <- .sl_at_level(hits, "Molecular subspecies", db)
  fallback  <- FALSE

  if (length(species) == 0) {
    siblings <- .get_many(substr(key, 1, 14), db$by_key14)
    fb <- .sl_at_level(siblings, "Species", db)
    if (length(fb) == 1) {
      species  <- fb
      fallback <- TRUE
    }
  }

  empty$sl_tier      <- as.integer(tier)
  empty$sl_n_hits    <- length(hits)
  empty$sl_species   <- if (length(species) == 0) NA_character_ else paste(species, collapse = "; ")
  empty$sl_molecular <- if (length(molecular) == 0) NA_character_ else paste(molecular, collapse = "; ")
  empty$sl_n_species <- length(species)
  empty$sl_fallback  <- fallback
  empty
}
