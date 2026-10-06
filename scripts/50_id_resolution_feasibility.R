# Feasibility test: can PubChem CIDs and HMDB IDs be turned into the
# right InChIKeys?
#
# Uses the development files (LEOCC Metabolon and Broad 2022), where rows
# carry both identifiers and an InChIKey, so the key stored in the file
# is the reference. The COMETS data are not touched here.
#
# Criteria are in notes/decisions.md (2026-10-06), written before this
# script was first run:
#   CID-derived keys agree with the stored key on the first 14
#   characters in at least 95% of comparable rows.
#   HMDB-derived keys: at least 90%.
#
# CIDs: all rows with exactly one CID and a valid stored key.
# HMDB: a random sample of 300 such rows (one request per ID).
#
# Needs: httr2, jsonlite, internet. Answers are cached, so a rerun does
# not ask PubChem again.
# Output: results/id_resolution/
# Run from the project root: Rscript scripts/50_id_resolution_feasibility.R

source("scripts/50_id_functions.R")

out_dir <- "results/id_resolution"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache_path <- file.path(out_dir, "pubchem_cache.rds")
cache <- list(env = if (file.exists(cache_path)) readRDS(cache_path) else new.env())
on.exit(saveRDS(cache$env, cache_path), add = TRUE)

base <- "data-raw/MetLinkR/inst/extdata/Example_output"
leocc <- read.csv(file.path(base, "LEOCC_Metabolon_Annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))
broad <- read.csv(file.path(base, "Broad_2022Aug_annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))

rows <- rbind(
  data.frame(row_id = paste0("LEOCC_Metabolon:", seq_len(nrow(leocc))),
             name = leocc$CHEMICAL_NAME, stored_key = leocc$INCHIKEY,
             cid_raw = leocc$PUBCHEM, hmdb_raw = leocc$HMDB),
  data.frame(row_id = paste0("Broad_2022Aug:", seq_len(nrow(broad))),
             name = broad$Metabolite, stored_key = broad$inChIKey,
             cid_raw = broad$pubChemId, hmdb_raw = broad$hmdbId)
)
valid_key <- !is.na(rows$stored_key) & grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", rows$stored_key)
rows$cid  <- single_cid(rows$cid_raw)
rows$hmdb <- single_hmdb(rows$hmdb_raw)
rows$cid_list  <- lengths(id_tokens(rows$cid_raw)) > 1
rows$hmdb_list <- lengths(id_tokens(rows$hmdb_raw)) > 1


# ---- CID route -----------------------------------------------------------

cid_rows <- rows[valid_key & !is.na(rows$cid), ]
cid_keys <- cid_to_key(cid_rows$cid, cache)
cid_rows$derived_key <- unname(cid_keys[cid_rows$cid])
saveRDS(cache$env, cache_path)


# ---- HMDB route ----------------------------------------------------------
# A HMDB ID can map to several CIDs. It counts as resolved only when all
# of them give the same InChIKey.

hmdb_pool <- rows[valid_key & !is.na(rows$hmdb), ]
set.seed(20261006)
hmdb_rows <- hmdb_pool[sample(nrow(hmdb_pool), min(300, nrow(hmdb_pool))), ]
hmdb_cids <- lapply(hmdb_rows$hmdb, hmdb_to_cids, cache = cache)
saveRDS(cache$env, cache_path)
all_cids <- unique(unlist(hmdb_cids))
hmdb_cid_keys <- cid_to_key(all_cids, cache)
saveRDS(cache$env, cache_path)
hmdb_rows$n_cids <- lengths(hmdb_cids)
hmdb_rows$derived_key <- vapply(hmdb_cids, function(cs) {
  k <- unique(stats::na.omit(unname(hmdb_cid_keys[cs])))
  if (length(k) == 1) k else NA_character_
}, character(1))
hmdb_rows$status <- ifelse(hmdb_rows$n_cids == 0, "no_cid",
                    ifelse(is.na(hmdb_rows$derived_key), "several_keys", "resolved"))


# ---- Compare -------------------------------------------------------------

compare <- function(d) {
  d$agree_14 <- !is.na(d$derived_key) & substr(d$derived_key, 1, 14) == substr(d$stored_key, 1, 14)
  d$agree_25 <- !is.na(d$derived_key) & substr(d$derived_key, 1, 25) == substr(d$stored_key, 1, 25)
  d
}
cid_rows  <- compare(cid_rows)
hmdb_rows <- compare(hmdb_rows)
write.csv(cid_rows,  file.path(out_dir, "cid_route.csv"),  row.names = FALSE)
write.csv(hmdb_rows, file.path(out_dir, "hmdb_route.csv"), row.names = FALSE)

pct <- function(a, b) round(100 * a / b, 1)
c_res <- sum(!is.na(cid_rows$derived_key))
h_res <- sum(!is.na(hmdb_rows$derived_key))

message("Rows skipped because the cell lists several IDs: CID ", sum(rows$cid_list & valid_key),
        ", HMDB ", sum(rows$hmdb_list & valid_key))
message("\nCID route: ", nrow(cid_rows), " rows, ", c_res, " resolved")
message("  agree on first 14: ", sum(cid_rows$agree_14), " of ", c_res,
        " (", pct(sum(cid_rows$agree_14), c_res), "%, criterion 95%)")
message("  agree on first 25: ", sum(cid_rows$agree_25), " of ", c_res)
message("\nHMDB route: ", nrow(hmdb_rows), " sampled rows")
print(table(hmdb_rows$status))
message("  agree on first 14: ", sum(hmdb_rows$agree_14), " of ", h_res,
        " (", pct(sum(hmdb_rows$agree_14), h_res), "%, criterion 90%)")
message("  agree on first 25: ", sum(hmdb_rows$agree_25), " of ", h_res)

show <- c("row_id", "name", "stored_key", "derived_key")
message("\nCID rows that disagree on the first 14 characters:")
print(cid_rows[!is.na(cid_rows$derived_key) & !cid_rows$agree_14, c(show, "cid")], row.names = FALSE)
message("\nHMDB rows that disagree on the first 14 characters:")
print(hmdb_rows[!is.na(hmdb_rows$derived_key) & !hmdb_rows$agree_14, c(show, "hmdb")], row.names = FALSE)
message("\nFiles written to ", out_dir)
