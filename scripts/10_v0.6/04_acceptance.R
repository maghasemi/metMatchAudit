# Acceptance checks for v0.6: LMSD as a second lipid source and the
# Biocrates p180 vendor table (notes/decisions.md, 2026-10-07).
#
# Development data: every pair whose rule changed against v0.4.2 is
# listed, the evaluation against notes/review.csv is recomputed, and new
# lipid match proposals are written out for review by hand.
# COMETS: after adaptation. Changed pairs are listed with the expert's
# view, and new lipid matches are counted against the expert groups.
#
# Run from the project root: Rscript scripts/10_v0.6/04_acceptance.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.6.0")

out_dir <- "results/v06_check"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
options(width = 250)

lmsd_rds <- "data-raw/lipidmaps/lmsd_db.rds"
lmsd <- if (file.exists(lmsd_rds)) readRDS(lmsd_rds) else {
  db <- load_lmsd("data-raw/lipidmaps/structures.sdf", expected_md5 = "134c2b5c6d76e9883d46339a6289eefd")
  saveRDS(db, lmsd_rds)
  db
}
print(lmsd)
sl <- load_swisslipids("data-raw/swisslipids/lipids.tsv", expected_md5 = "b6ddd5bc3003967dce95e905e5748953")

rerun <- function(old, cache_from, label, files = old$info$files) {
  cache <- file.path(out_dir, paste0("pubchem_cache_", label, ".rds"))
  if (!file.exists(cache)) file.copy(cache_from, cache)
  new <- audit_metlinkr(files, swisslipids = sl, cache_file = cache, lmsd = lmsd)
  saveRDS(new, file.path(out_dir, paste0("audit_", label, ".rds")))
  new
}

changed_pairs <- function(old, new) {
  m <- merge(old$pairs[, c("row_a", "row_b", "group_id", "relationship", "rule_id")],
             new$pairs[, c("row_a", "row_b", "relationship", "rule_id", "evidence")],
             by = c("row_a", "row_b"), suffixes = c("_old", "_new"))
  ch <- m[m$rule_id_old != m$rule_id_new, ]
  nm <- stats::setNames(new$rows$original_name, new$rows$row_id)
  ch$name_a <- nm[ch$row_a]; ch$name_b <- nm[ch$row_b]
  ch
}


# ---- Development data ------------------------------------------------------------

old_dev <- readRDS("results/dev_v0.4.2/audit.rds")  # v0.4.2 baseline, see decisions.md 2026-10-07
dev <- rerun(old_dev, "results/pubchem_cache.rds", "development")

message("\n==== development")
message("Lipid sources:"); print(table(dev$rows$lipid_source, useNA = "ifany"))
ch <- changed_pairs(old_dev, dev)
message("Pairs whose rule changed against v0.4.2: ", nrow(ch))
if (nrow(ch) > 0) print(ch[, c("name_a", "name_b", "rule_id_old", "rule_id_new", "relationship_new")],
                       row.names = FALSE, right = FALSE)
write.csv(ch, file.path(out_dir, "dev_changed_pairs.csv"), row.names = FALSE)

review <- read.csv("notes/review.csv", colClasses = "character")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")
g <- dev$groups[match(review$metlinkr_name, dev$groups$group_id), ]
rel <- g$worst_relationship
outcome <- ifelse(review$category == "E", ifelse(rel == "exact", "E_exact", "E_not_exact"),
           ifelse(rel == unname(expected_map[review$category]), "agree",
           ifelse(rel == "unresolved", "safe_miss", "wrong")))
message("Evaluation against the review (v0.4.2: 35 agree, 8 safe_miss, 0 wrong, 10 E_not_exact):")
print(table(outcome))

old_lm <- paste(old_dev$lipid_matches$row_a, old_dev$lipid_matches$row_b)
new_lm <- dev$lipid_matches[!paste(dev$lipid_matches$row_a, dev$lipid_matches$row_b) %in% old_lm, ]
message("Lipid match proposals: ", nrow(old_dev$lipid_matches), " before, ", nrow(dev$lipid_matches),
        " now, ", nrow(new_lm), " new")
new_lm$correct <- ""; new_lm$note <- ""
write.csv(new_lm, file.path(out_dir, "dev_new_lipid_matches_for_review.csv"), row.names = FALSE)


# ---- COMETS, after adaptation -----------------------------------------------------

old_c <- readRDS("results/comets_v0.4.2/audit.rds")
files <- old_c$info$files
files$vendor <- ifelse(files$label == "inputs_1539353822", "biocrates_p180", NA_character_)
com <- rerun(old_c, "results/comets_v0.4.2/pubchem_cache.rds", "comets", files)

# Expert groups, joined as in scripts/08_comets/02_blind_test.R.
src <- "data-raw/comets"
norm_id <- function(x) toupper(gsub("[[:space:]]", "", x))
rows <- com$rows
ids <- unlist(lapply(files$label, function(s)
  read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)[[1]]))
rk <- paste(rows$file, norm_id(ids), sep = "|")
mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character", na.strings = c("", "NA"))
mk <- paste(sub("\\.xlsx$", "", mm$SOURCE_FILE), norm_id(mm$IPT_METABID), sep = "|")
groups_of <- split(mm$Manual_Metabolite_Group_Label, mk)[rk]
names(groups_of) <- rows$row_id
same <- function(a, b) vapply(seq_along(a), function(i)
  length(intersect(groups_of[[a[i]]], groups_of[[b[i]]])) > 0, logical(1))

message("\n==== comets (after adaptation)")
message("Lipid sources:"); print(table(rows$lipid_source, useNA = "ifany"))
ch <- changed_pairs(old_c, com)
ch$expert_same <- same(ch$row_a, ch$row_b)
message("Pairs whose rule changed against v0.4.2: ", nrow(ch))
if (nrow(ch) > 0) print(ch[order(ch$relationship_new), c("name_a", "name_b", "rule_id_old", "rule_id_new",
                                                         "relationship_new", "expert_same")],
                       row.names = FALSE, right = FALSE)
write.csv(ch, file.path(out_dir, "comets_changed_pairs.csv"), row.names = FALSE)

lm <- com$lipid_matches
lm <- lm[!lm$already_grouped, ]
lm$expert_same <- same(lm$row_a, lm$row_b)
message("New lipid matches across files, not grouped by MetLinkR: ", nrow(lm),
        ", in the same expert group: ", sum(lm$expert_same))
print(table(lm$relationship, lm$expert_same))
write.csv(lm, file.path(out_dir, "comets_new_lipid_matches.csv"), row.names = FALSE)
