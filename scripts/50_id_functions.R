# Identifier parsing and PubChem lookups for turning PubChem CIDs and
# HMDB IDs into InChIKeys. Used by scripts/50_id_resolution_feasibility.R
# and meant to move into the package once the test passes.

# All tokens in an identifier cell, with prefixes such as "pubchem:" or
# "hmdb:" removed. Cells can hold lists ("hmdb:HMDB0002189,HMDB0001276").
id_tokens <- function(x) {
  lapply(as.character(x), function(cell) {
    if (is.na(cell) || !nzchar(trimws(cell))) return(character(0))
    parts <- trimws(unlist(strsplit(cell, "[,;| ]+")))
    parts <- sub("^(pubchem|hmdb|cid)[:_]", "", parts, ignore.case = TRUE)
    parts[nzchar(parts)]
  })
}

# One PubChem CID per cell, or NA when the cell is empty, holds several
# IDs, or holds something that is not a number.
single_cid <- function(x) {
  vapply(id_tokens(x), function(t) {
    t <- unique(t)
    if (length(t) == 1 && grepl("^[0-9]+$", t)) t else NA_character_
  }, character(1))
}

# One HMDB ID per cell in the current seven-digit form. "HMDB00641" and
# "HMDB641" both become "HMDB0000641". NA for empty cells and lists.
single_hmdb <- function(x) {
  vapply(id_tokens(x), function(t) {
    t <- unique(toupper(t))
    if (length(t) != 1 || !grepl("^HMDB[0-9]+$", t)) return(NA_character_)
    digits <- sub("^HMDB", "", t)
    if (nchar(digits) > 7) return(NA_character_)
    paste0("HMDB", formatC(as.integer(digits), width = 7, flag = "0"))
  }, character(1))
}

# Cached GET or POST to PUG REST. Returns the parsed JSON, or a list with
# an error message.
pug <- function(path, body = NULL, cache) {
  key <- paste(path, if (is.null(body)) "" else paste(unlist(body), collapse = ","))
  if (!is.null(cache$env[[key]])) return(cache$env[[key]])
  req <- httr2::request(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/", path)) |>
    httr2::req_user_agent("metMatchAudit") |>
    httr2::req_timeout(60) |>
    httr2::req_error(is_error = function(r) FALSE)
  if (!is.null(body)) req <- httr2::req_body_form(req, !!!body)
  resp <- tryCatch(httr2::req_perform(req), error = function(e) e)
  out <- if (inherits(resp, "error")) {
    list(error = conditionMessage(resp))
  } else {
    parsed <- tryCatch(jsonlite::fromJSON(httr2::resp_body_string(resp)),
                       error = function(e) list(error = "not JSON"))
    parsed$http_status <- httr2::resp_status(resp)
    parsed
  }
  cache$env[[key]] <- out
  Sys.sleep(0.25)  # stay under PubChem's 5 requests per second
  out
}

# InChIKey for each CID, in batches. CIDs that PubChem does not return
# stay NA.
cid_to_key <- function(cids, cache, batch = 100) {
  cids <- unique(cids[!is.na(cids)])
  out <- stats::setNames(rep(NA_character_, length(cids)), cids)
  for (start in seq(1, length(cids), by = batch)) {
    chunk <- cids[start:min(start + batch - 1, length(cids))]
    res <- pug("compound/cid/property/InChIKey/JSON",
               body = list(cid = paste(chunk, collapse = ",")), cache = cache)
    props <- res$PropertyTable$Properties
    if (!is.null(props) && nrow(props) > 0) {
      out[as.character(props$CID)] <- props$InChIKey
    }
  }
  out
}

# CIDs PubChem lists for an HMDB registry ID.
hmdb_to_cids <- function(hmdb, cache) {
  res <- pug(paste0("compound/xref/RegistryID/", hmdb, "/cids/JSON"), cache = cache)
  cids <- res$IdentifierList$CID
  if (is.null(cids)) character(0) else as.character(cids)
}
