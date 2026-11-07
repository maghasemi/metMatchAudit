# Acceptance check for v0.9.1 (notes/decisions.md, 2026-11-07, v0.9.1 plan).
#
# Reruns the development data (with the Broad name fallback) and COMETS on
# copies of the v0.9 PubChem caches. When no request fails, the pairs must
# be identical to results/v09_check/, and no request may fail.
#
# Every answer should already be in the copied caches, so this should need
# no internet. If a cache held a failed answer from an older version, it is
# asked again, and that request shows in the count.
#
# Run from the project root:
#   Rscript scripts/17_v0.9.1/01_acceptance.R 2>&1 | tee results/v091_check_log.txt

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.9.1")
options(width = 200)

out <- "results/v091_check"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

lmsd <- readRDS("data-raw/lipidmaps/lmsd_db.rds")
stopifnot(identical(lmsd$md5, "134c2b5c6d76e9883d46339a6289eefd"))
sl <- load_swisslipids("data-raw/swisslipids/lipids.tsv", expected_md5 = "b6ddd5bc3003967dce95e905e5748953")

rerun <- function(files, cache_from, cache_to) {
  if (!file.exists(cache_to)) file.copy(cache_from, cache_to)
  audit_metlinkr(files, swisslipids = sl, cache_file = cache_to, lmsd = lmsd)
}

pair_cols <- c("group_id", "row_a", "row_b", "relationship", "rule_id", "evidence")
same_pairs <- function(old, new) identical(as.list(old$pairs[, pair_cols]), as.list(new$pairs[, pair_cols]))
same_keys <- function(old, new) identical(old$rows[, c("row_id", "inchikey", "inchikey_source")],
                                          new$rows[, c("row_id", "inchikey", "inchikey_source")])

old_dev <- readRDS("results/v09_check/audit_development.rds")
old_com <- readRDS("results/v09_check/audit_comets_no_fallback.rds")
stopifnot(old_dev$info$package_version == "0.9.0", old_com$info$package_version == "0.9.0")

dev <- rerun(old_dev$info$files, "results/v09_check/pubchem_cache_development.rds",
             file.path(out, "pubchem_cache_development.rds"))
com <- rerun(old_com$info$files, "results/v09_check/pubchem_cache_comets.rds",
             file.path(out, "pubchem_cache_comets.rds"))
saveRDS(dev, file.path(out, "audit_development.rds"))
saveRDS(com, file.path(out, "audit_comets.rds"))

summary <- data.frame(
  check = c("development: failed PubChem requests", "development: keys identical to 0.9.0",
            "development: pairs identical to 0.9.0",
            "COMETS: failed PubChem requests", "COMETS: keys identical to 0.9.0",
            "COMETS: pairs identical to 0.9.0"),
  value = c(dev$info$pubchem_failed, same_keys(old_dev, dev), same_pairs(old_dev, dev),
            com$info$pubchem_failed, same_keys(old_com, com), same_pairs(old_com, com)))
summary$pass <- c(dev$info$pubchem_failed == 0, summary$value[2:3] == 1,
                  com$info$pubchem_failed == 0, summary$value[5:6] == 1)
message("Development pairs: ", nrow(old_dev$pairs), " in 0.9.0, ", nrow(dev$pairs), " now")
message("COMETS pairs: ", nrow(old_com$pairs), " in 0.9.0, ", nrow(com$pairs), " now")
message("\n==== summary")
print(summary, row.names = FALSE, right = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
message("\nAcceptance: ", if (isTRUE(all(summary$pass))) "ALL PASS" else "AT LEAST ONE FAILED")
