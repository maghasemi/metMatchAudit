# Third blind test, step 6: runs the audit once on A, C1 and C2 and
# compares it with the committed labels.
#
# Plan and criteria: notes/decisions.md (2026-10-07, third blind test).
# Label md5s and the stereo rule: notes/blind3_labels.md. Both were
# committed before this script was first run.
#
# Privacy: file A is private. Everything that names an A row goes to
# data-raw/own/blind3/results/, which is gitignored. Only counts go to
# results/blind3/, which is committed.
#
# Run from the project root: Rscript scripts/13_blind_test_3/04_test.R

pkgload::load_all(".", quiet = TRUE)
options(width = 200)

# ---- Frozen code and fixed labels -----------------------------------------

stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.7.0")
if (system("git diff --quiet v0.7.0 -- R DESCRIPTION NAMESPACE inst", ignore.stdout = TRUE) != 0) {
  stop("Package code differs from tag v0.7.0. The test must run on the frozen version.")
}
note <- "notes/blind3_labels.md"
if (system(paste("git ls-files --error-unmatch", note), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0 ||
    system(paste("git diff --quiet HEAD --", note)) != 0) {
  stop(note, " must be committed, unchanged, before the run.")
}

lab  <- "data-raw/own/blind3"
priv <- file.path(lab, "results")
pub  <- "results/blind3"
if (file.exists(file.path(priv, "audit.rds"))) {
  stop("The audit has run already. The test runs once; delete ", priv, " only on purpose.")
}
dir.create(priv, showWarnings = FALSE, recursive = TRUE)
dir.create(pub,  showWarnings = FALSE, recursive = TRUE)

md5 <- function(p) unname(tools::md5sum(p))
expected <- c(label_pairs.csv  = "898c093671f12638e61cf3637661ec3d",
              label_lipids.csv = "1a769432ed8b2c9d54d5f461dc2eeca2")
for (f in names(expected)) {
  got <- md5(file.path(lab, f))
  if (!identical(got, expected[[f]])) stop(f, " has md5 ", got, ", expected ", expected[[f]], ".")
}
stopifnot(md5("data-raw/own/A_metabolon.csv") == "db6a4afbbbcd2cdbe5015a5f234c154e",
          md5("data-raw/lai2021/C1_lai_feces.csv") == "8bb3413252fb77f36db75dc88b187632",
          md5("data-raw/lai2021/C2_lai_serum.csv") == "1e8400cbfbcb13f7440901e8426fa037")

lp <- read.csv(file.path(lab, "label_pairs.csv"),  colClasses = "character", na.strings = "")
ll <- read.csv(file.path(lab, "label_lipids.csv"), colClasses = "character", na.strings = "")
stopifnot(nrow(lp) == 94, all(lp$category %in% LETTERS[1:6]), nrow(ll) == 40)


# ---- The audit, once -----------------------------------------------------

mout <- "data-raw/own/metlinkr_output_blind3"
files <- data.frame(
  path         = file.path(mout, c("A_metabolon_metLinkR.csv", "C1_lai_feces_metLinkR.csv",
                                   "C2_lai_serum_metLinkR.csv")),
  name_col     = "BIOCHEMICAL",
  inchikey_col = NA,
  pubchem_col  = "PUBCHEM",
  hmdb_col     = c("HMDB_ID", NA, NA),
  label        = c("A", "C1", "C2"),
  stringsAsFactors = FALSE
)
lmsd <- readRDS("data-raw/lipidmaps/lmsd_db.rds")
stopifnot(identical(lmsd$md5, "134c2b5c6d76e9883d46339a6289eefd"))

audit <- audit_metlinkr(files,
                        swisslipids  = "data-raw/swisslipids/lipids.tsv",
                        expected_md5 = "b6ddd5bc3003967dce95e905e5748953",
                        cache_file   = file.path(priv, "pubchem_cache.rds"),
                        lmsd         = lmsd)
saveRDS(audit, file.path(priv, "audit.rds"))
print(audit)
stopifnot(sum(audit$rows$file == "A") == 494, sum(audit$rows$file == "C1") == 533,
          sum(audit$rows$file == "C2") == 231)


# ---- Scoring ---------------------------------------------------------------

key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))
lab_order <- c("E", "C", "F", "D", "B", "A")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")

i <- match(key(lp$row_a, lp$row_b), key(audit$pairs$row_a, audit$pairs$row_b))
if (anyNA(i)) stop("Some labeled pairs are missing from the audit.")
lp$tool <- audit$pairs$relationship[i]
lp$rule <- audit$pairs$rule_id[i]
lp$with_A <- grepl("^A:", lp$row_a) | grepl("^A:", lp$row_b)
lp$definite <- lp$tool %in% c("exact", "distinct")
lp$c1 <- ifelse(!lp$definite, NA,
         ifelse(lp$tool == "exact", lp$category == "A", lp$category %in% c("C", "E")))

group_score <- function(p) {
  gl <- tapply(p$category, p$group_id, function(x) lab_order[min(match(x, lab_order))])
  g <- data.frame(group_id = names(gl), category = unname(gl), stringsAsFactors = FALSE)
  j <- match(g$group_id, audit$groups$group_id)
  if (anyNA(j)) stop("Some labeled groups are missing from the audit.")
  g$tool <- audit$groups$worst_relationship[j]
  g$rule <- audit$groups$worst_rule[j]
  g$outcome <- with(g, ifelse(category == "E", ifelse(tool == "exact", "E_exact", "E_not_exact"),
                       ifelse(tool == expected_map[category], "agree",
                       ifelse(tool == "unresolved", "safe_miss", "wrong"))))
  g
}
groups <- group_score(lp)
in_c2 <- groups$category %in% c("A", "B", "C", "D")

# C3: proposals across files that involve a sampled row; truth pairs from
# the aligned ";" lists in the lipid sheet.
truth <- do.call(rbind, lapply(seq_len(nrow(ll)), function(k) {
  m <- trimws(strsplit(ll$matches_in_other_files[k], ";")[[1]])
  m <- m[nzchar(m) & m != "none"]
  if (length(m) == 0) return(NULL)
  data.frame(row = ll$row[k], other = m, stringsAsFactors = FALSE)
}))
truth_keys <- unique(key(truth$row, truth$other))
lm <- audit$lipid_matches
lm <- lm[lm$file_a != lm$file_b & (lm$row_a %in% ll$row | lm$row_b %in% ll$row), ]
lm$in_labels <- key(lm$row_a, lm$row_b) %in% truth_keys
found <- truth_keys %in% key(lm$row_a, lm$row_b)

pct <- function(x, n) if (n == 0) NA_real_ else round(100 * x / n, 1)
c1_right <- sum(lp$c1, na.rm = TRUE); c1_total <- sum(lp$definite)
c2_wrong <- sum(groups$outcome[in_c2] == "wrong"); c2_total <- sum(in_c2)
summary <- data.frame(
  criterion = c("C1: definite claims that agree with the label (at least 90%)",
                "C2: groups labeled A to D with a wrong definite answer (at most 5%)",
                "C3: lipid proposals that are pairs in the labels (at least 90%)",
                "C3 recall: labeled lipid pairs the tool proposed (no threshold)",
                "C1 on pairs without A (reported, no threshold)"),
  count = c(c1_right, c2_wrong, sum(lm$in_labels), sum(found),
            sum(lp$c1[!lp$with_A], na.rm = TRUE)),
  of    = c(c1_total, c2_total, nrow(lm), length(truth_keys), sum(lp$definite[!lp$with_A])),
  stringsAsFactors = FALSE)
summary$percent <- mapply(pct, summary$count, summary$of)
summary$met <- c(if (c1_total == 0) NA else c1_right / c1_total >= 0.9,
                 c2_wrong / c2_total <= 0.05,
                 if (nrow(lm) == 0) NA else mean(lm$in_labels) >= 0.9,
                 NA, NA)

message("\n==== Third blind test")
print(summary, row.names = FALSE, right = FALSE)
message("\nPairs, label against tool:")
print(table(label = lp$category, tool = factor(lp$tool, levels = .severity)))
message("\nPairs, label and tool by rule:")
print(table(paste(lp$category, lp$tool), lp$rule))
message("\nGroups, label against outcome:")
print(table(label = groups$category, outcome = groups$outcome))
message("\nLipid proposals involving sampled rows, by tool relationship:")
print(table(tool = lm$relationship, in_labels = lm$in_labels))

# Counts only: committed.
write.csv(summary, file.path(pub, "summary.csv"), row.names = FALSE)
write.csv(subset(as.data.frame(table(label = lp$category, tool = lp$tool, rule = lp$rule, with_A = lp$with_A)), Freq > 0),
          file.path(pub, "pairs_by_label_tool_rule.csv"), row.names = FALSE)
write.csv(subset(as.data.frame(table(label = groups$category, outcome = groups$outcome)), Freq > 0),
          file.path(pub, "groups_by_label_outcome.csv"), row.names = FALSE)
writeLines(c(paste("package:", audit$info$package_version),
             paste("run:", audit$info$run_date),
             paste("swisslipids md5:", audit$info$swisslipids_md5),
             paste("lmsd md5:", audit$info$lmsd_md5),
             paste("structure checks:", audit$info$structure_checks),
             paste("metlinkr output md5:", basename(files$path), md5(files$path)),
             paste("label md5:", names(expected), expected)),
           file.path(pub, "run_info.txt"))

# Names and rows: private.
write.csv(lp, file.path(priv, "pairs_compared.csv"), row.names = FALSE, na = "")
write.csv(groups, file.path(priv, "groups_compared.csv"), row.names = FALSE, na = "")
write.csv(lm, file.path(priv, "lipid_proposals_compared.csv"), row.names = FALSE, na = "")
write.csv(data.frame(key = truth_keys, found = found), file.path(priv, "lipid_labels_found.csv"), row.names = FALSE)
write_audit(audit, file.path(priv, "audit"))
message("\nCounts in ", pub, "; details in ", priv, " (not committed).")
