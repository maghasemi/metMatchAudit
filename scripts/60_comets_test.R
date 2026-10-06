# COMETS external test of metMatchAudit v0.4.0.
#
# Compares MetLinkR alone, and MetLinkR plus the audit, with the expert
# mappings of the COMETS consortium (ManualMappings.csv, MetLinkR SI) on
# five datasets that played no part in developing the tool.
#
# Criteria and definitions are in notes/decisions.md (2026-10-06, COMETS
# test), committed before this script was first run. Run it once.
#
# Needs: the MetLinkR output from metlinkr-run/run_metlinkr_comets.R in
# data-raw/comets/metlinkr_output/, the SwissLipids file, internet for
# identifier resolution on the first run.
# Output: results/comets_test/
# Run from the project root: Rscript scripts/60_comets_test.R

pkgload::load_all(".", quiet = TRUE)

# ---- The package must be the frozen v0.4.0 --------------------------------
# Commits after the tag (criteria, this script) are fine, but nothing in
# the package code may differ from it.
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.4.0")
changed <- system("git diff --quiet v0.4.0 -- R DESCRIPTION NAMESPACE", ignore.stdout = TRUE)
if (changed != 0) stop("Package code differs from tag v0.4.0. The test must run on the frozen version.")

out_dir <- "results/comets_test"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

src  <- "data-raw/comets"
mout <- file.path(src, "metlinkr_output")
stems <- c("inputs_1539353822", "inputs_1539354673", "inputs_1553799316",
           "inputs_1555609921", "inputs_1571898296")
nightingale <- "inputs_1571898296"


# ---- Run the audit on MetLinkR's output -----------------------------------

files <- data.frame(
  path         = file.path(mout, paste0(stems, "_metLinkR.csv")),
  name_col     = "metabolite_name",
  inchikey_col = NA,
  pubchem_col  = c(NA, NA, NA, "PUBCHEM", NA),
  hmdb_col     = c(NA, "HMDB_ID", NA, NA, NA),
  label        = stems,
  stringsAsFactors = FALSE
)
audit <- audit_metlinkr(files,
                        swisslipids  = "data-raw/swisslipids/lipids.tsv",
                        expected_md5 = "b6ddd5bc3003967dce95e905e5748953",
                        cache_file   = file.path(out_dir, "pubchem_cache.rds"))
saveRDS(audit, file.path(out_dir, "audit.rds"))
rows <- audit$rows


# ---- Attach the expert group to every row ---------------------------------
# Rows are matched by file and metabolite ID. The match must be one to one;
# anything else stops the test, since a partial match would bias it.

mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character",
               na.strings = c("", "NA"))
mm$stem <- sub("\\.xlsx$", "", mm$SOURCE_FILE)
mm_key <- paste(mm$stem, mm$IPT_METABID, sep = "|")
if (anyDuplicated(mm_key)) stop("ManualMappings.csv has duplicate file and ID pairs.")

row_metabid <- unlist(lapply(stems, function(s) {
  d <- read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)
  d[[1]]  # first column is the metabolite ID in all five files (metabid or METABID)
}))
if (length(row_metabid) != nrow(rows)) stop("Input files and audit rows differ in number.")
rows$expert_group <- mm$Manual_Metabolite_Group_Label[match(paste(rows$file, row_metabid, sep = "|"), mm_key)]

unmatched <- sum(is.na(rows$expert_group))
if (unmatched > 0) {
  print(table(rows$file[is.na(rows$expert_group)]))
  stop(unmatched, " rows have no expert group. Check the ID columns before going on.")
}
if (nrow(mm) != nrow(rows)) stop("ManualMappings.csv and the input files differ in number of rows.")


# ---- Helpers -------------------------------------------------------------------

# Number of pairs of rows from different files within each group.
cross_pairs_in <- function(group, file) {
  ok <- !is.na(group)
  n_g  <- table(group[ok])
  n_gf <- table(paste(group[ok], file[ok]))
  sum(choose(n_g, 2)) - sum(choose(n_gf, 2))
}
expert_of <- stats::setNames(rows$expert_group, rows$row_id)
file_of   <- stats::setNames(rows$file, rows$row_id)

prf <- function(tp, predicted, positives) {
  p <- if (predicted == 0) NA else tp / predicted
  r <- tp / positives
  f <- if (is.na(p) || (p + r) == 0) NA else 2 * p * r / (p + r)
  c(true_positive = tp, predicted = predicted, precision = round(p, 3),
    recall = round(r, 3), f1 = round(f, 3))
}


# ---- Pairs ----------------------------------------------------------------------

positives <- cross_pairs_in(rows$expert_group, rows$file)

pairs <- audit$pairs
pairs$cross_file <- file_of[pairs$row_a] != file_of[pairs$row_b]
pairs <- pairs[pairs$cross_file, ]
pairs$expert_same <- expert_of[pairs$row_a] == expert_of[pairs$row_b]
pairs$nightingale <- file_of[pairs$row_a] == nightingale | file_of[pairs$row_b] == nightingale
write.csv(pairs, file.path(out_dir, "cross_file_pairs.csv"), row.names = FALSE)

# Question 1: are the tool's definite claims right?
exact    <- pairs[pairs$relationship == "exact", ]
distinct <- pairs[pairs$relationship == "distinct", ]
definite_right <- sum(exact$expert_same) + sum(!distinct$expert_same)
definite_total <- nrow(exact) + nrow(distinct)
by_relationship <- as.data.frame.matrix(table(pairs$relationship,
                                              ifelse(pairs$expert_same, "expert_same", "expert_different")))

# Question 2: MetLinkR alone against MetLinkR plus the audit.
metlinkr_alone <- prf(sum(pairs$expert_same), nrow(pairs), positives)
kept <- pairs[pairs$relationship != "distinct", ]
with_audit <- prf(sum(kept$expert_same), nrow(kept), positives)
strict <- pairs[pairs$relationship %in% c("exact", "broader_narrower"), ]
with_audit_strict <- prf(sum(strict$expert_same), nrow(strict), positives)
comparison <- rbind(metlinkr_alone = metlinkr_alone, with_audit = with_audit,
                    with_audit_exact_or_bn_only = with_audit_strict)

# Question 3: lipid matches MetLinkR did not make.
lm <- audit$lipid_matches
lm <- lm[!lm$already_grouped, ]
lm$expert_same <- expert_of[lm$row_a] == expert_of[lm$row_b]
write.csv(lm, file.path(out_dir, "new_lipid_matches.csv"), row.names = FALSE)

# Group level, as reported in the MetLinkR paper: share of expert groups
# spanning two or more files whose rows all got one MetLinkR name.
eg <- split(rows$group_id, rows$expert_group)
ef <- split(rows$file, rows$expert_group)
cross_groups <- names(ef)[vapply(ef, function(f) length(unique(f)) > 1, logical(1))]
found <- vapply(eg[cross_groups], function(g) !anyNA(g) && length(unique(g)) == 1, logical(1))

# Nightingale apart.
ng <- pairs[pairs$nightingale, ]

# Every wrong definite claim, for review afterwards.
nm <- stats::setNames(rows$original_name, rows$row_id)
wrong <- rbind(exact[!exact$expert_same, ], distinct[distinct$expert_same, ])
wrong$name_a <- nm[wrong$row_a]; wrong$name_b <- nm[wrong$row_b]
write.csv(wrong[, c("group_id", "name_a", "name_b", "relationship", "rule_id", "evidence")],
          file.path(out_dir, "wrong_definite_claims.csv"), row.names = FALSE)


# ---- Criteria (decisions.md, 2026-10-06) -----------------------------------------

c1 <- definite_right / definite_total
c2 <- with_audit["precision"] >= metlinkr_alone["precision"] &&
      isTRUE(with_audit["f1"] >= metlinkr_alone["f1"])
summary_table <- data.frame(
  criterion = c("C1: definite claims agree with the expert (at least 90%)",
                "C2: with the audit, precision not lower and F1 not lower"),
  value = c(sprintf("%d of %d (%.1f%%)", definite_right, definite_total, 100 * c1),
            sprintf("precision %.3f vs %.3f, F1 %.3f vs %.3f",
                    with_audit["precision"], metlinkr_alone["precision"],
                    with_audit["f1"], metlinkr_alone["f1"])),
  met = c(c1 >= 0.90, c2)
)
write.csv(summary_table, file.path(out_dir, "criteria.csv"), row.names = FALSE)
write.csv(comparison, file.path(out_dir, "comparison.csv"))
write.csv(by_relationship, file.path(out_dir, "by_relationship.csv"))

message("Rows: ", nrow(rows), "; expert cross-file positive pairs: ", positives)
message("Cross-file pairs inside MetLinkR groups: ", nrow(pairs))
message("\nCriteria:"); print(summary_table, row.names = FALSE)
message("\nTool relationship against the expert:"); print(by_relationship)
message("\nMetLinkR alone and with the audit:"); print(comparison)
message("\nNew lipid matches: ", nrow(lm), ", in the same expert group: ", sum(lm$expert_same))
message("Expert cross-file groups found whole by MetLinkR: ", sum(found), " of ", length(found),
        sprintf(" (%.1f%%; the paper reports 85.3%%)", 100 * mean(found)))
message("Nightingale pairs: ", nrow(ng), ", expert same: ", sum(ng$expert_same))
message("Key sources:"); print(audit$info$keys_by_source)
message("\nFiles written to ", out_dir)
