# Retest of identifier resolution with the revised design, on data not
# used so far: annotations.csv from the MetLinkR SI. That file is not
# one of the five files with expert mappings, and only identifier
# resolution is tested here, not any grouping.
#
# Rows are kept only if neither their CID nor the first 14 characters of
# their key occur in the development files (LEOCC, Broad), since the
# revisions were designed after seeing disagreements there.
#
# Criteria (notes/decisions.md, 2026-10-06, retest), same thresholds as
# the first run:
#   CID route: at least 95% agree with the stored key on the first 14.
#   HMDB route (random sample of 300, seed 20261006): at least 90%.
#
# Needs: httr2, jsonlite, internet; the package itself (load_all) for the
# name checks. Output: results/id_resolution_retest/
# Run from the project root: Rscript scripts/51_id_resolution_retest.R

pkgload::load_all(".", quiet = TRUE)
source("scripts/50_id_functions.R")

out_dir <- "results/id_resolution_retest"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache_path <- file.path(out_dir, "pubchem_cache.rds")
cache <- list(env = if (file.exists(cache_path)) readRDS(cache_path) else new.env())

ann <- read.csv("data-raw/comets/annotations.csv", colClasses = "character",
                na.strings = c("", "NA"))
base <- "data-raw/MetLinkR/inst/extdata/Example_output"
leocc <- read.csv(file.path(base, "LEOCC_Metabolon_Annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))
broad <- read.csv(file.path(base, "Broad_2022Aug_annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))
dev_cids <- unique(c(unlist(id_tokens(leocc$PUBCHEM)), unlist(id_tokens(broad$pubChemId))))
dev_b1   <- unique(substr(c(leocc$INCHIKEY, broad$inChIKey), 1, 14))

ann$row_id <- paste0("annotations:", seq_len(nrow(ann)))
valid_key <- !is.na(ann$inChIKey) & grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", ann$inChIKey)
ann$cid  <- single_cid(ann$pubChemId)
ann$hmdb <- single_hmdb(ann$hmdbId)
fresh <- valid_key & !(ann$cid %in% dev_cids) & !(substr(ann$inChIKey, 1, 14) %in% dev_b1)


# ---- CID route -----------------------------------------------------------

cid_rows <- ann[fresh & !is.na(ann$cid), c("row_id", "Metabolite", "inChIKey", "cid")]
res <- resolve_cid_keys(cid_rows$cid, cid_rows$Metabolite, cache)
saveRDS(cache$env, cache_path)
cid_rows$derived_key <- res$key
cid_rows$skipped_name <- res$skipped_name


# ---- HMDB route ----------------------------------------------------------

hmdb_pool <- ann[fresh & !is.na(ann$hmdb), c("row_id", "Metabolite", "inChIKey", "hmdb")]
set.seed(20261006)
hmdb_rows <- hmdb_pool[sample(nrow(hmdb_pool), min(300, nrow(hmdb_pool))), ]
skip_h <- is_group_name(hmdb_rows$Metabolite) | has_alternative_names(hmdb_rows$Metabolite)
hmdb_cids <- lapply(seq_len(nrow(hmdb_rows)), function(i) {
  if (skip_h[i]) character(0) else hmdb_to_cids(hmdb_rows$hmdb[i], cache)
})
saveRDS(cache$env, cache_path)
parents <- cid_to_parent(unique(unlist(hmdb_cids)), cache)
pkeys   <- cid_to_key(unique(unname(parents)), cache)
saveRDS(cache$env, cache_path)
hmdb_rows$skipped_name <- skip_h
hmdb_rows$derived_key <- vapply(hmdb_cids, function(cs) {
  if (length(cs) == 0) return(NA_character_)
  k <- unique(stats::na.omit(unname(pkeys[parents[cs]])))
  if (length(k) == 1) k else NA_character_
}, character(1))


# ---- Compare -------------------------------------------------------------

agree <- function(d) {
  d$agree_14 <- !is.na(d$derived_key) & substr(d$derived_key, 1, 14) == substr(d$inChIKey, 1, 14)
  d$agree_25 <- !is.na(d$derived_key) & substr(d$derived_key, 1, 25) == substr(d$inChIKey, 1, 25)
  d
}
cid_rows <- agree(cid_rows); hmdb_rows <- agree(hmdb_rows)
write.csv(cid_rows,  file.path(out_dir, "cid_route.csv"),  row.names = FALSE)
write.csv(hmdb_rows, file.path(out_dir, "hmdb_route.csv"), row.names = FALSE)

pct <- function(a, b) round(100 * a / b, 1)
c_res <- sum(!is.na(cid_rows$derived_key)); h_res <- sum(!is.na(hmdb_rows$derived_key))

message("Fresh rows (CID and key block not in the development files): ", sum(fresh))
message("\nCID route: ", nrow(cid_rows), " rows, ", sum(cid_rows$skipped_name),
        " skipped for group or alternative names, ", c_res, " resolved")
message("  agree on first 14: ", sum(cid_rows$agree_14), " of ", c_res,
        " (", pct(sum(cid_rows$agree_14), c_res), "%, criterion 95%)")
message("  agree on first 25: ", sum(cid_rows$agree_25), " of ", c_res)
message("\nHMDB route: ", nrow(hmdb_rows), " sampled, ", sum(hmdb_rows$skipped_name),
        " skipped for names, ", h_res, " resolved")
message("  agree on first 14: ", sum(hmdb_rows$agree_14), " of ", h_res,
        " (", pct(sum(hmdb_rows$agree_14), h_res), "%, criterion 90%)")

show <- c("row_id", "Metabolite", "inChIKey", "derived_key")
message("\nCID rows that disagree on the first 14 characters:")
print(cid_rows[!is.na(cid_rows$derived_key) & !cid_rows$agree_14, c(show, "cid")], row.names = FALSE)
message("\nHMDB rows that disagree on the first 14 characters:")
print(hmdb_rows[!is.na(hmdb_rows$derived_key) & !hmdb_rows$agree_14, c(show, "hmdb")], row.names = FALSE)
message("\nFiles written to ", out_dir)
