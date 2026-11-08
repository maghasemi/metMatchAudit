# Acceptance checks for v0.11 (notes/decisions.md, 2026-11-07, v0.11 plan).
# All results are after adaptation.
#
# Reruns the development data (with the Broad name fallback) and COMETS on
# copies of the v0.9.1 PubChem caches and compares with results/v091_check/:
#   1. the only pairs that change are R01 pairs, and all move to R01s;
#   2. development evaluation: 35 agree, 8 safe_miss, 0 wrong, 0 E_exact,
#      10 E_not_exact, with same_name counted as a safe miss;
#   3. no lipid match proposal is lost or gained;
#   4. COMETS definite claims (exact, distinct) across files are unchanged;
#   5. COMETS same_name pairs the expert covers: at least 98% in one expert
#      group, and every disagreement listed.
# The unit tests are run separately, before this script.
#
# Run from the project root:
#   Rscript scripts/18_v0.11/02_acceptance.R 2>&1 | tee results/v11_check_log.txt

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.11.0")
options(width = 200)

out <- "results/v11_check"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

lmsd <- readRDS("data-raw/lipidmaps/lmsd_db.rds")
stopifnot(identical(lmsd$md5, "134c2b5c6d76e9883d46339a6289eefd"))
sl <- load_swisslipids("data-raw/swisslipids/lipids.tsv", expected_md5 = "b6ddd5bc3003967dce95e905e5748953")

rerun <- function(files, cache_from, cache_to) {
  if (!file.exists(cache_to)) file.copy(cache_from, cache_to)
  audit_metlinkr(files, swisslipids = sl, cache_file = cache_to, lmsd = lmsd)
}

old_dev <- readRDS("results/v091_check/audit_development.rds")
old_com <- readRDS("results/v091_check/audit_comets.rds")
stopifnot(old_dev$info$package_version == "0.9.1", old_com$info$package_version == "0.9.1")

dev <- rerun(old_dev$info$files, "results/v091_check/pubchem_cache_development.rds",
             file.path(out, "pubchem_cache_development.rds"))
com <- rerun(old_com$info$files, "results/v091_check/pubchem_cache_comets.rds",
             file.path(out, "pubchem_cache_comets.rds"))
saveRDS(dev, file.path(out, "audit_development.rds"))
saveRDS(com, file.path(out, "audit_comets.rds"))
message("Failed PubChem requests: development ", dev$info$pubchem_failed, ", COMETS ", com$info$pubchem_failed)


# ---- 1. Which pairs changed --------------------------------------------------

changed <- function(old, new) {
  m <- merge(old$pairs[, c("row_a", "row_b", "group_id", "relationship", "rule_id")],
             new$pairs[, c("row_a", "row_b", "relationship", "rule_id")],
             by = c("row_a", "row_b"), suffixes = c("_old", "_new"))
  stopifnot(nrow(m) == nrow(old$pairs), nrow(m) == nrow(new$pairs))
  m[m$rule_id_old != m$rule_id_new | m$relationship_old != m$relationship_new, ]
}
ch_dev <- changed(old_dev, dev)
ch_com <- changed(old_com, com)
only_r01 <- function(ch) all(ch$rule_id_old == "R01" & ch$rule_id_new == "R01s" & ch$relationship_new == "same_name")
message("\n==== 1. changed pairs")
message("Development: ", nrow(ch_dev), "; COMETS: ", nrow(ch_com))
print(table(set = c(rep("development", nrow(ch_dev)), rep("comets", nrow(ch_com))),
            change = c(paste(ch_dev$rule_id_old, "->", ch_dev$rule_id_new),
                       paste(ch_com$rule_id_old, "->", ch_com$rule_id_new))))
ok_changes <- only_r01(ch_dev) && only_r01(ch_com)
message("Only R01 pairs changed, all to R01s: ", ok_changes)


# ---- 2. Development evaluation ------------------------------------------------

review <- read.csv("notes/review.csv", colClasses = "character")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")
rel <- dev$groups$worst_relationship[match(review$metlinkr_name, dev$groups$group_id)]
outcome <- ifelse(review$category == "E", ifelse(rel == "exact", "E_exact", "E_not_exact"),
           ifelse(rel == unname(expected_map[review$category]), "agree",
           ifelse(rel %in% c("unresolved", "same_name"), "safe_miss", "wrong")))
counts <- table(factor(outcome, levels = c("agree", "safe_miss", "wrong", "E_exact", "E_not_exact")))
message("\n==== 2. development evaluation (0.9.1: 35 agree, 8 safe_miss, 0 wrong, 0 E_exact, 10 E_not_exact)")
print(counts)
ok_eval <- counts[["agree"]] == 35 && counts[["safe_miss"]] == 8 && counts[["wrong"]] == 0 &&
  counts[["E_exact"]] == 0 && counts[["E_not_exact"]] == 10


# ---- 3. Lipid match proposals ------------------------------------------------

lm_key <- function(a) paste(a$lipid_matches$row_a, a$lipid_matches$row_b, a$lipid_matches$relationship)
lm_same <- function(old, new) setequal(lm_key(old), lm_key(new))
message("\n==== 3. lipid match proposals")
message("Development: ", nrow(old_dev$lipid_matches), " before, ", nrow(dev$lipid_matches), " now; same set: ",
        lm_same(old_dev, dev))
message("COMETS: ", nrow(old_com$lipid_matches), " before, ", nrow(com$lipid_matches), " now; same set: ",
        lm_same(old_com, com))
ok_lipids <- lm_same(old_dev, dev) && lm_same(old_com, com)


# ---- 4. COMETS definite claims across files ----------------------------------

file_of <- stats::setNames(com$rows$file, com$rows$row_id)
definite <- function(a) {
  p <- a$pairs[a$pairs$relationship %in% c("exact", "distinct"), ]
  p <- p[file_of[p$row_a] != file_of[p$row_b], ]
  sort(paste(p$row_a, p$row_b, p$relationship))
}
ok_definite <- identical(definite(old_com), definite(com))
message("\n==== 4. COMETS definite claims across files: ", length(definite(old_com)), " before, ",
        length(definite(com)), " now; identical: ", ok_definite)


# ---- 5. COMETS same_name pairs against the expert ----------------------------

# The same join as scripts/12_v0.7/02_acceptance.R. Empty labels are
# dropped, and a row without a label counts as not covered.
src <- "data-raw/comets"
norm_id <- function(x) toupper(gsub("[[:space:]]", "", x))
ids <- unlist(lapply(com$info$files$label, function(s)
  read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)[[1]]))
stopifnot(length(ids) == nrow(com$rows))
rk <- paste(com$rows$file, norm_id(ids), sep = "|")
mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character", na.strings = c("", "NA"))
mm <- mm[!is.na(mm$Manual_Metabolite_Group_Label), ]
mk <- paste(sub("\\.xlsx$", "", mm$SOURCE_FILE), norm_id(mm$IPT_METABID), sep = "|")
groups_of <- split(mm$Manual_Metabolite_Group_Label, mk)[rk]
names(groups_of) <- com$rows$row_id

sn <- com$pairs[com$pairs$relationship == "same_name", c("row_a", "row_b", "group_id", "evidence")]
sn$covered <- lengths(groups_of[sn$row_a]) > 0 & lengths(groups_of[sn$row_b]) > 0
sn$expert_same <- vapply(seq_len(nrow(sn)), function(i)
  length(intersect(groups_of[[sn$row_a[i]]], groups_of[[sn$row_b[i]]])) > 0, logical(1))
known <- sn[sn$covered, ]
share <- if (nrow(known) > 0) mean(known$expert_same) else NA_real_
message("\n==== 5. COMETS same_name pairs: ", nrow(sn), "; the expert covers ", nrow(known),
        "; in one expert group: ", sum(known$expert_same), " (", round(100 * share, 1), "%). Needed: at least 98%.")
nm <- stats::setNames(com$rows$original_name, com$rows$row_id)
if (any(!known$expert_same)) {
  message("same_name pairs the expert separates:")
  d <- known[!known$expert_same, ]
  print(data.frame(name_a = nm[d$row_a], name_b = nm[d$row_b], group_id = d$group_id), row.names = FALSE, right = FALSE)
}
sn$name_a <- unname(nm[sn$row_a]); sn$name_b <- unname(nm[sn$row_b])
write.csv(sn, file.path(out, "comets_same_name_pairs.csv"), row.names = FALSE)
write.csv(rbind(cbind(set = rep("development", nrow(ch_dev)), ch_dev), cbind(set = rep("comets", nrow(ch_com)), ch_com)),
          file.path(out, "changed_pairs.csv"), row.names = FALSE)


# ---- Summary -----------------------------------------------------------------

summary <- data.frame(
  check = c("no failed PubChem request", "only R01 pairs change, all to R01s", "development evaluation as in 0.9.1",
            "no lipid proposal lost or gained", "COMETS definite claims unchanged",
            "COMETS same_name pairs in one expert group"),
  value = c(dev$info$pubchem_failed + com$info$pubchem_failed, ok_changes, ok_eval, ok_lipids, ok_definite,
            round(share, 4)),
  pass = c(dev$info$pubchem_failed + com$info$pubchem_failed == 0, ok_changes, ok_eval, ok_lipids, ok_definite,
           isTRUE(share >= 0.98)))
message("\n==== summary")
print(summary, row.names = FALSE, right = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
message("\nAutomatic checks: ", if (isTRUE(all(summary$pass))) "ALL PASS" else "AT LEAST ONE FAILED")
