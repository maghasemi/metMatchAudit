# Acceptance checks for v0.9 (notes/decisions.md, 2026-11-07, v0.9 plan and
# the feasibility entry after it). All results are after adaptation.
#
#   1. Without name_fallback_col: the pairs of the development data and of
#      COMETS must be identical to the saved runs in results/v07_check/.
#      Those runs were made with 0.7.0; 0.7.1 changed no behaviour.
#   2. With name_fallback_col = "name" on Broad: evaluation against
#      notes/review.csv, every pair whose rule changed, lipid match
#      proposals lost and gained.
#   3. A sheet for the hand check of every changed pair and every new
#      lipid match proposal.
#
# Every run gets a copy of the PubChem cache of the run it is compared
# with, so only the package differs. Only public example data and COMETS
# counts are written here. No private file is read.
#
# Run from the project root:
#   Rscript scripts/15_v0.9/02_acceptance.R 2>&1 | tee results/v09_check_log.txt

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.9.0")
options(width = 220)

out <- "results/v09_check"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

lmsd <- readRDS("data-raw/lipidmaps/lmsd_db.rds")
stopifnot(identical(lmsd$md5, "134c2b5c6d76e9883d46339a6289eefd"))
sl <- load_swisslipids("data-raw/swisslipids/lipids.tsv", expected_md5 = "b6ddd5bc3003967dce95e905e5748953")

rerun <- function(files, cache_from, cache_to) {
  if (!file.exists(cache_to)) file.copy(cache_from, cache_to)
  audit_metlinkr(files, swisslipids = sl, cache_file = cache_to, lmsd = lmsd)
}

pair_cols <- c("group_id", "row_a", "row_b", "relationship", "rule_id", "evidence")
same_pairs <- function(old, new) {
  identical(as.list(old$pairs[, pair_cols]), as.list(new$pairs[, pair_cols]))
}

changed_pairs <- function(old, new) {
  m <- merge(old$pairs[, c("row_a", "row_b", "group_id", "relationship", "rule_id")],
             new$pairs[, c("row_a", "row_b", "relationship", "rule_id", "evidence")],
             by = c("row_a", "row_b"), suffixes = c("_old", "_new"))
  ch <- m[m$rule_id_old != m$rule_id_new, ]
  nm <- stats::setNames(new$rows$original_name, new$rows$row_id)
  ch$name_a <- unname(nm[ch$row_a]); ch$name_b <- unname(nm[ch$row_b])
  ch
}

lipid_diff <- function(old, new) {
  ko <- paste(old$lipid_matches$row_a, old$lipid_matches$row_b)
  kn <- paste(new$lipid_matches$row_a, new$lipid_matches$row_b)
  list(lost = old$lipid_matches[!ko %in% kn, ], gained = new$lipid_matches[!kn %in% ko, ])
}

show <- function(x, cols) {
  if (nrow(x) == 0) return(invisible())
  print(x[, cols], row.names = FALSE, right = FALSE)
}


# ---- 1. Without name_fallback_col --------------------------------------------

old_dev <- readRDS("results/v07_check/audit_development.rds")
old_com <- readRDS("results/v07_check/audit_comets.rds")
cache_dev <- file.path(out, "pubchem_cache_development.rds")
cache_com <- file.path(out, "pubchem_cache_comets.rds")

dev0 <- rerun(old_dev$info$files, "results/v07_check/pubchem_cache_development.rds", cache_dev)
com0 <- rerun(old_com$info$files, "results/v07_check/pubchem_cache_comets.rds", cache_com)
saveRDS(dev0, file.path(out, "audit_development_no_fallback.rds"))
saveRDS(com0, file.path(out, "audit_comets_no_fallback.rds"))

same_dev <- same_pairs(old_dev, dev0)
same_com <- same_pairs(old_com, com0)
message("==== 1. without name_fallback_col")
message("Development pairs: ", nrow(old_dev$pairs), " before, ", nrow(dev0$pairs),
        " now; identical: ", same_dev)
message("COMETS pairs: ", nrow(old_com$pairs), " before, ", nrow(com0$pairs),
        " now; identical: ", same_com)
if (!same_dev) show(changed_pairs(old_dev, dev0), c("name_a", "name_b", "rule_id_old", "rule_id_new"))
if (!same_com) show(changed_pairs(old_com, com0), c("name_a", "name_b", "rule_id_old", "rule_id_new"))


# ---- 2. With name_fallback_col on Broad --------------------------------------

files <- old_dev$info$files
files$name_fallback_col <- ifelse(files$label == "Broad_2022Aug", "name", NA_character_)
dev <- rerun(files, "results/v07_check/pubchem_cache_development.rds", cache_dev)
saveRDS(dev, file.path(out, "audit_development.rds"))

message("\n==== 2. with name_fallback_col = \"name\" on Broad")
message("Rows by name source:")
print(table(file = dev$rows$file, name_source = dev$rows$name_source, useNA = "ifany"))

filled <- dev$rows$row_id[dev$rows$name_source %in% "fallback"]
fr <- dev$rows[dev$rows$row_id %in% filled, ]
message("Filled rows: ", nrow(fr), "; with a lipid profile: ", sum(!is.na(fr$lipid_species)),
        "; with a valid key: ", sum(fr$ik_valid))
print(table(level = fr$lipid_level, source = fr$lipid_source, useNA = "ifany"))
if (any(is.na(fr$lipid_species))) {
  message("Filled rows without a lipid profile:")
  show(fr[is.na(fr$lipid_species), ], c("row_id", "original_name"))
}

message("\nPairs by relationship, before and after:")
rel_levels <- c("exact", "broader_narrower", "stereo_variant", "distinct", "unresolved")
print(rbind(before = table(factor(dev0$pairs$relationship, levels = rel_levels)),
            after  = table(factor(dev$pairs$relationship, levels = rel_levels))))

ch <- changed_pairs(dev0, dev)
only_filled <- all(ch$row_a %in% filled | ch$row_b %in% filled)
same_pair_set <- nrow(dev0$pairs) == nrow(dev$pairs)
message("\nPairs whose rule changed: ", nrow(ch), ". Expected from the feasibility reading: ",
        "10 to R05 and 3 to R08.")
print(table(change = paste(ch$rule_id_old, "->", ch$rule_id_new)))
show(ch, c("name_a", "name_b", "rule_id_old", "rule_id_new", "relationship_new"))
message("All changed pairs hold a filled row: ", only_filled)
message("Same set of pairs before and after: ", same_pair_set)
write.csv(ch, file.path(out, "dev_changed_pairs.csv"), row.names = FALSE)

review <- read.csv("notes/review.csv", colClasses = "character")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")
rel <- dev$groups$worst_relationship[match(review$metlinkr_name, dev$groups$group_id)]
outcome <- ifelse(review$category == "E", ifelse(rel == "exact", "E_exact", "E_not_exact"),
           ifelse(rel == unname(expected_map[review$category]), "agree",
           ifelse(rel == "unresolved", "safe_miss", "wrong")))
outcome <- factor(outcome, levels = c("agree", "safe_miss", "wrong", "E_exact", "E_not_exact"))
counts <- table(outcome)
message("\nEvaluation against the review (v0.7: 35 agree, 8 safe_miss, 0 wrong, 0 E_exact, 10 E_not_exact):")
print(counts)

ld <- lipid_diff(dev0, dev)
message("\nLipid match proposals: ", nrow(dev0$lipid_matches), " before, ", nrow(dev$lipid_matches),
        " now; lost ", nrow(ld$lost), ", gained ", nrow(ld$gained))
lm_cols <- c("name_a", "name_b", "lipid_species", "relationship", "rule_id", "shape", "already_grouped")
show(ld$lost, lm_cols)
show(ld$gained, lm_cols)
write.csv(ld$lost, file.path(out, "dev_lipid_matches_lost.csv"), row.names = FALSE)
write.csv(ld$gained, file.path(out, "dev_lipid_matches_gained.csv"), row.names = FALSE)

message("\nRow flags by check, before and after:")
flag_levels <- sort(unique(c(dev0$row_flags$check, dev$row_flags$check)))
print(rbind(before = table(factor(dev0$row_flags$check, levels = flag_levels)),
            after  = table(factor(dev$row_flags$check, levels = flag_levels))))


# ---- 3. Sheet for the hand check ---------------------------------------------

hand <- rbind(
  data.frame(kind = rep("changed_pair", nrow(ch)), row_a = ch$row_a, row_b = ch$row_b,
             name_a = ch$name_a, name_b = ch$name_b, rule_old = ch$rule_id_old,
             rule = ch$rule_id_new, relationship = ch$relationship_new, evidence = ch$evidence,
             stringsAsFactors = FALSE),
  data.frame(kind = rep("lipid_proposal", nrow(ld$gained)), row_a = ld$gained$row_a, row_b = ld$gained$row_b,
             name_a = ld$gained$name_a, name_b = ld$gained$name_b, rule_old = rep(NA_character_, nrow(ld$gained)),
             rule = ld$gained$rule_id, relationship = ld$gained$relationship, evidence = ld$gained$evidence,
             stringsAsFactors = FALSE))
hand$right <- ""
hand$note <- ""
sheet <- file.path(out, "hand_check.csv")
message("\n==== 3. hand check")
if (file.exists(sheet)) {
  message("Left as it is, so nothing you filled in is lost: ", sheet)
} else {
  write.csv(hand, sheet, row.names = FALSE, na = "")
  message("Wrote ", nrow(hand), " lines to ", sheet, ". Fill `right` with yes or no on every line.")
}


# ---- Summary -----------------------------------------------------------------

summary <- data.frame(
  check = c("development pairs identical without fallback", "COMETS pairs identical without fallback",
            "wrong stays at 0", "no source-error group exact", "agree at least 35",
            "only pairs with a filled row change rule", "no lipid proposal lost"),
  value = c(same_dev, same_com, counts[["wrong"]], counts[["E_exact"]], counts[["agree"]],
            only_filled, nrow(ld$lost)),
  pass = c(same_dev, same_com, counts[["wrong"]] == 0, counts[["E_exact"]] == 0, counts[["agree"]] >= 35,
           only_filled, nrow(ld$lost) == 0))
message("\n==== summary")
print(summary, row.names = FALSE, right = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
message("\nAutomatic checks: ", if (all(summary$pass)) "ALL PASS" else "AT LEAST ONE FAILED",
        ". The hand check is still open: ", nrow(hand), " lines, at least 95% must be right.")
