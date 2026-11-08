# Fourth blind test, step 6: runs the audit once on F1, F2 and F3 and
# compares it with the committed labels.
#
# Plan and criteria C1 to C4: notes/decisions.md (2026-11-08, fourth blind
# test, plan). Label md5s and the stereo rule: notes/blind4_labels.md. Both
# were committed before this script was first run.
#
# The files are CC BY-NC-ND. Everything that names a row goes to
# data-raw/blind4/results/, which is gitignored. Only counts go to
# results/blind4/, which is committed.
#
# Run from the project root:
#   Rscript scripts/19_blind_test_4/04_test.R 2>&1 | tee results/blind4_log.txt

pkgload::load_all(".", quiet = TRUE)
options(width = 200)

# ---- Frozen code and fixed labels -----------------------------------------

stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.11.0")
if (system("git diff --quiet v0.11.0 -- R DESCRIPTION NAMESPACE inst", ignore.stdout = TRUE) != 0) {
  stop("Package code differs from tag v0.11.0. The test must run on the frozen version.")
}
note <- "notes/blind4_labels.md"
if (system(paste("git ls-files --error-unmatch", note), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0 ||
    system(paste("git diff --quiet HEAD --", note)) != 0) {
  stop(note, " must be committed, unchanged, before the run.")
}

lab  <- "data-raw/blind4/labels"
priv <- "data-raw/blind4/results"
pub  <- "results/blind4"
if (file.exists(file.path(priv, "audit.rds"))) {
  stop("The audit has run already. The test runs once; delete ", priv, " only on purpose.")
}
dir.create(priv, showWarnings = FALSE, recursive = TRUE)
dir.create(pub,  showWarnings = FALSE, recursive = TRUE)

md5 <- function(p) unname(tools::md5sum(p))
expected <- c(label_pairs.csv  = "f9f03def8a5b04546c8839673ea6d811",
              label_lipids.csv = "70e2213390b427e7d58adf24ec81a051")
for (f in names(expected)) {
  got <- md5(file.path(lab, f))
  if (!identical(got, expected[[f]])) stop(f, " has md5 ", got, ", expected ", expected[[f]], ".")
}
mout <- "data-raw/blind4/metlinkr_output"
stopifnot(md5("data-raw/blind4/F1.csv") == "4f1009d43f2512718999809135555ee6",
          md5("data-raw/blind4/F2.csv") == "fe888f7912fb44e3182f9ec91427784d",
          md5("data-raw/blind4/F3.csv") == "6d7c70b4a7ccd256ece6b28361944131",
          md5(file.path(mout, "F1_metLinkR.csv")) == "fd0905bdee002465fbd5b0c9ef86e482",
          md5(file.path(mout, "F2_metLinkR.csv")) == "b1924621969128019f5f1895044d0fc1",
          md5(file.path(mout, "F3_metLinkR.csv")) == "4157b66b6012d9bf95d051dd8f75f506")

lp <- read.csv(file.path(lab, "label_pairs.csv"),  colClasses = "character", na.strings = "")
ll <- read.csv(file.path(lab, "label_lipids.csv"), colClasses = "character", na.strings = "")
stopifnot(nrow(lp) == 67, all(lp$category %in% LETTERS[1:6]), nrow(ll) == 40,
          all(lp$stratum %in% c("same", "different")))


# ---- The audit, once -----------------------------------------------------

files <- data.frame(
  path         = file.path(mout, c("F1_metLinkR.csv", "F2_metLinkR.csv", "F3_metLinkR.csv")),
  name_col     = "BIOCHEMICAL",
  inchikey_col = NA,
  pubchem_col  = c("PUBCHEM", NA, NA),
  hmdb_col     = c("HMDB", NA, NA),
  label        = c("F1", "F2", "F3"),
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
stopifnot(sum(audit$rows$file == "F1") == 991, sum(audit$rows$file == "F2") == 919,
          sum(audit$rows$file == "F3") == 463)


# ---- Scoring ---------------------------------------------------------------

key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))
lab_order <- c("E", "C", "F", "D", "B", "A")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")

i <- match(key(lp$row_a, lp$row_b), key(audit$pairs$row_a, audit$pairs$row_b))
if (anyNA(i)) stop("Some labeled pairs are missing from the audit.")
lp$tool <- audit$pairs$relationship[i]
lp$rule <- audit$pairs$rule_id[i]
lp$definite <- lp$tool %in% c("exact", "distinct")
lp$c1 <- ifelse(!lp$definite, NA,
         ifelse(lp$tool == "exact", lp$category == "A", lp$category %in% c("C", "E")))

# same_name counts as unresolved: a safe miss, never a wrong answer.
group_score <- function(p) {
  gl <- tapply(p$category, p$group_id, function(x) lab_order[min(match(x, lab_order))])
  g <- data.frame(group_id = names(gl), category = unname(gl), stringsAsFactors = FALSE)
  g$stratum <- p$stratum[match(g$group_id, p$group_id)]
  j <- match(g$group_id, audit$groups$group_id)
  if (anyNA(j)) stop("Some labeled groups are missing from the audit.")
  g$tool <- audit$groups$worst_relationship[j]
  g$rule <- audit$groups$worst_rule[j]
  g$outcome <- with(g, ifelse(category == "E", ifelse(tool == "exact", "E_exact", "E_not_exact"),
                       ifelse(tool == expected_map[category], "agree",
                       ifelse(tool %in% c("unresolved", "same_name"), "safe_miss", "wrong"))))
  g
}
groups <- group_score(lp)
groups$in_c2 <- groups$category %in% c("A", "B", "C", "D")

# C3: proposals across files that involve a sampled row. Truth pairs come
# from the ";" lists in the lipid sheet (A and D matches both count).
truth <- do.call(rbind, lapply(seq_len(nrow(ll)), function(k) {
  m <- trimws(strsplit(ll$matches_in_other_files[k], ";")[[1]])
  m <- m[!is.na(m) & nzchar(m) & m != "none"]
  if (length(m) == 0) return(NULL)
  data.frame(row = ll$row[k], other = m, stringsAsFactors = FALSE)
}))
truth_keys <- unique(key(truth$row, truth$other))
lm <- audit$lipid_matches
lm <- lm[lm$file_a != lm$file_b & (lm$row_a %in% ll$row | lm$row_b %in% ll$row), ]
lm$in_labels <- key(lm$row_a, lm$row_b) %in% truth_keys
found <- truth_keys %in% key(lm$row_a, lm$row_b)

pct <- function(x, n) if (n == 0) NA_real_ else round(100 * x / n, 1)
score <- function(p, g, label) {
  c1_total <- sum(p$definite); c1_right <- sum(p$c1, na.rm = TRUE)
  g2 <- g[g$in_c2, ]
  c2_wrong <- sum(g2$outcome == "wrong"); c2_total <- nrow(g2)
  sn <- p$tool == "same_name"; sn_a <- sum(sn & p$category == "A")
  data.frame(
    part = label,
    c1_right = c1_right, c1_of = c1_total, c1_pct = pct(c1_right, c1_total),
    c2_wrong = c2_wrong, c2_of = c2_total, c2_pct = pct(c2_wrong, c2_total),
    c4_a = sn_a, c4_of = sum(sn), c4_pct = pct(sn_a, sum(sn)),
    stringsAsFactors = FALSE)
}
by_part <- rbind(score(lp, groups, "all"),
                 score(lp[lp$stratum == "same", ], groups[groups$stratum == "same", ], "same-name stratum"),
                 score(lp[lp$stratum == "different", ], groups[groups$stratum == "different", ], "different-names stratum"))
by_part$c1_met <- ifelse(by_part$c1_of == 0, NA, by_part$c1_right / by_part$c1_of >= 0.9)
by_part$c2_met <- ifelse(by_part$c2_of == 0, NA, by_part$c2_wrong / by_part$c2_of <= 0.05)
by_part$c4_met <- ifelse(by_part$c4_of == 0, NA, by_part$c4_a / by_part$c4_of >= 0.95)
by_part$c4_note <- ifelse(by_part$c4_of < 20, "fewer than 20 pairs: too few to carry the threshold", "")

c3 <- data.frame(
  part = c("C3: lipid proposals that are pairs in the labels (at least 90%)",
           "C3 recall: labeled lipid pairs the tool proposed (no threshold)"),
  count = c(sum(lm$in_labels), sum(found)),
  of = c(nrow(lm), length(truth_keys)), stringsAsFactors = FALSE)
c3$percent <- mapply(pct, c3$count, c3$of)
c3$met <- c(if (nrow(lm) == 0) NA else mean(lm$in_labels) >= 0.9, NA)

message("\n==== Fourth blind test: C1, C2, C4, all pairs and by stratum")
print(by_part, row.names = FALSE, right = FALSE)
message("\n==== C3")
print(c3, row.names = FALSE, right = FALSE)
message("\nPairs, label against tool:")
print(table(label = lp$category, tool = factor(lp$tool, levels = .severity)))
message("\nPairs, label and tool by rule:")
print(table(paste(lp$category, lp$tool), lp$rule))
message("\nGroups, label against outcome:")
print(table(label = groups$category, outcome = groups$outcome))
message("\nLipid proposals involving sampled rows, by tool relationship:")
print(table(tool = lm$relationship, in_labels = lm$in_labels))
message("\nSampled lipid rows by file: ", paste(names(table(substr(ll$row, 1, 2))), table(substr(ll$row, 1, 2)), collapse = ", "))

# Counts only: committed.
write.csv(by_part, file.path(pub, "summary_c1_c2_c4.csv"), row.names = FALSE)
write.csv(c3, file.path(pub, "summary_c3.csv"), row.names = FALSE)
write.csv(subset(as.data.frame(table(stratum = lp$stratum, label = lp$category, tool = lp$tool, rule = lp$rule)), Freq > 0),
          file.path(pub, "pairs_by_label_tool_rule.csv"), row.names = FALSE)
write.csv(subset(as.data.frame(table(stratum = groups$stratum, label = groups$category, outcome = groups$outcome)), Freq > 0),
          file.path(pub, "groups_by_label_outcome.csv"), row.names = FALSE)
writeLines(c(paste("package:", audit$info$package_version),
             paste("run:", audit$info$run_date),
             paste("swisslipids md5:", audit$info$swisslipids_md5),
             paste("lmsd md5:", audit$info$lmsd_md5),
             paste("structure checks:", audit$info$structure_checks),
             paste("failed PubChem requests:", audit$info$pubchem_failed),
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
