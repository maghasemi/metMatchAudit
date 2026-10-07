# Acceptance checks for v0.5 (notes/decisions.md, 2026-10-07, v0.5 plan).
#
# Runs the audit with the RefMet check on the development data and on
# COMETS, with the same input files as the v0.4.2 runs, and checks:
#   - pair results are identical to v0.4.2, since no pair rule changed;
#   - the six rows from reviewed C and E groups get C06;
#   - every C06 flag is written to a sheet for review by hand.
# COMETS results are after adaptation and only counted.
#
# Run from the project root: Rscript scripts/71_v05_check.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.5.0")

refmet <- load_refmet("data-raw/refmet/refmet_all.txt",
                      expected_md5 = "63079d6af6aff5fec9ddc9d3f952afa2")
out_dir <- "results/v05_check"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

run <- function(old_path, cache_from, label) {
  old <- readRDS(old_path)
  cache <- file.path(out_dir, paste0("pubchem_cache_", label, ".rds"))
  if (!file.exists(cache)) file.copy(cache_from, cache)
  new <- audit_metlinkr(old$info$files,
                        swisslipids  = "data-raw/swisslipids/lipids.tsv",
                        expected_md5 = "b6ddd5bc3003967dce95e905e5748953",
                        cache_file   = cache, refmet = refmet)
  saveRDS(new, file.path(out_dir, paste0("audit_", label, ".rds")))

  cols <- c("row_a", "row_b", "relationship", "rule_id", "evidence")
  same <- nrow(old$pairs) == nrow(new$pairs) &&
    all(vapply(cols, function(k) identical(as.character(old$pairs[[k]]),
                                            as.character(new$pairs[[k]])), logical(1)))
  message("\n==== ", label)
  message("Pair results identical to v0.4.2: ", same)
  message("RefMet check per row:")
  print(table(new$rows$rm_check, useNA = "ifany"))
  message("By basis, for rows that were compared:")
  cmp <- new$rows$rm_check %in% c("agrees", "conflict", "identifier_conflict", "unchecked")
  print(table(new$rows$rm_basis[cmp], new$rows$rm_check[cmp]))

  r <- new$rows[new$rows$rm_check %in% c("conflict", "identifier_conflict"),
                c("row_id", "file", "original_name", "rm_entry", "rm_check", "rm_basis",
                  "inchikey_source", "rm_detail")]
  r <- r[order(r$rm_check, r$original_name), ]
  list(audit = new, flagged = r, pairs_same = same)
}

dev    <- run("results/audit_v0.1_sealed.rds", "results/pubchem_cache.rds", "development")
comets <- run("results/comets_v0.4.2/audit.rds", "results/comets_v0.4.2/pubchem_cache.rds", "comets")

# Review sheet for the development data: every C06 and C07, with an empty
# column for the hand review. Only C06 counts toward the 80% check.
sheet <- dev$flagged
sheet$real <- ""
sheet$note <- ""
write.csv(sheet, file.path(out_dir, "dev_flags_for_review.csv"), row.names = FALSE)
write.csv(comets$flagged, file.path(out_dir, "comets_flags_after_adaptation.csv"), row.names = FALSE)

required <- c("andro steroid monosulfate C19H28O6S (1)*", "(15:0)-anacardic acid",
              "docosapentaenoate (n3 DPA; 22:5n3)", "N-acetyl-1-methylhistidine*",
              "octadecadienedioate (C18:2-DC)*", "pregnenediol sulfate (C21H34O5S)*")
got <- dev$audit$rows[dev$audit$rows$original_name %in% required &
                        dev$audit$rows$file == "LEOCC_Metabolon",
                      c("original_name", "rm_entry", "rm_check", "rm_basis")]
message("\nRequired rows (all must be conflict):")
print(got, row.names = FALSE, right = FALSE)
message("Required rows with C06: ", sum(got$rm_check == "conflict"), " of ", length(required))

options(width = 250)
message("\nDevelopment flags (also in ", file.path(out_dir, "dev_flags_for_review.csv"), "):")
print(dev$flagged[, c("original_name", "rm_entry", "rm_check", "rm_basis", "rm_detail")],
      row.names = FALSE, right = FALSE)
message("\nCOMETS flags, after adaptation:")
print(comets$flagged[, c("file", "original_name", "rm_entry", "rm_check", "rm_basis")],
      row.names = FALSE, right = FALSE)
