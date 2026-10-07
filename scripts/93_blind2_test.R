# Second blind test, step 6: runs the audit once on A and B and compares
# it with the committed labels.
#
# Plan and criteria: notes/decisions.md (2026-10-07, second blind test).
# Label md5s and the group label order: notes/blind2_labels.md. Both were
# committed before this script was first run.
#
# Privacy: file A is private. Everything that names an A row goes to
# data-raw/own/blind2/results/, which is gitignored. Only counts go to
# results/blind2/, which is committed.
#
# Needs: MetLinkR output in data-raw/own/metlinkr_output/, SwissLipids,
# LMSD, internet for identifier resolution.
# Run from the project root: Rscript scripts/93_blind2_test.R

pkgload::load_all(".", quiet = TRUE)
options(width = 200)

# ---- Frozen code and fixed labels -----------------------------------------

stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.6.0")
changed <- system("git diff --quiet v0.6.0 -- R DESCRIPTION NAMESPACE inst", ignore.stdout = TRUE)
if (changed != 0) stop("Package code differs from tag v0.6.0. The test must run on the frozen version.")

note <- "notes/blind2_labels.md"
if (system(paste("git ls-files --error-unmatch", note), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0 ||
    system(paste("git diff --quiet HEAD --", note)) != 0) {
  stop(note, " must be committed, unchanged, before the run.")
}

own  <- "data-raw/own"
lab  <- file.path(own, "blind2")
priv <- file.path(lab, "results")
pub  <- "results/blind2"
if (file.exists(file.path(priv, "audit.rds"))) {
  stop("The audit has run already. The test runs once; delete ", priv, " only on purpose.")
}
dir.create(priv, showWarnings = FALSE, recursive = TRUE)
dir.create(pub,  showWarnings = FALSE, recursive = TRUE)

md5 <- function(p) unname(tools::md5sum(p))
expected <- c(label_pairs.csv  = "3f70bedbd9a0c5a1f21204571a37a75d",
              label_lipids.csv = "bd170037bd8d528ddfc5d984fc5e9813")
for (f in names(expected)) {
  got <- md5(file.path(lab, f))
  if (!identical(got, expected[[f]])) stop(f, " has md5 ", got, ", expected ", expected[[f]], ".")
}

lp <- read.csv(file.path(lab, "label_pairs.csv"),  colClasses = "character", na.strings = "")
ll <- read.csv(file.path(lab, "label_lipids.csv"), colClasses = "character", na.strings = "")
stopifnot(nrow(lp) == 86, all(lp$category %in% LETTERS[1:6]),
          nrow(ll) == 40, all(!is.na(ll$matches_in_A)))


# ---- The audit, once -----------------------------------------------------

mout <- file.path(own, "metlinkr_output")
files <- data.frame(
  path         = file.path(mout, c("A_metabolon_metLinkR.csv", "B_ding2021_metLinkR.csv")),
  name_col     = "BIOCHEMICAL",
  inchikey_col = c(NA, "INCHIKEY"),
  smiles_col   = c(NA, "SMILES"),
  pubchem_col  = "PUBCHEM",
  hmdb_col     = c("HMDB_ID", NA),
  label        = c("A", "B"),
  stringsAsFactors = FALSE
)

lmsd_rds <- "data-raw/lipidmaps/lmsd_db.rds"
lmsd <- if (file.exists(lmsd_rds)) readRDS(lmsd_rds) else
  load_lmsd("data-raw/lipidmaps/structures.sdf", expected_md5 = "134c2b5c6d76e9883d46339a6289eefd")
stopifnot(identical(lmsd$md5, "134c2b5c6d76e9883d46339a6289eefd"))

audit <- audit_metlinkr(files,
                        swisslipids  = "data-raw/swisslipids/lipids.tsv",
                        expected_md5 = "b6ddd5bc3003967dce95e905e5748953",
                        cache_file   = file.path(priv, "pubchem_cache.rds"),
                        lmsd         = lmsd)
saveRDS(audit, file.path(priv, "audit.rds"))
print(audit)
stopifnot(sum(audit$rows$file == "A") == 494, sum(audit$rows$file == "B") == 1547)


# ---- C1, pairs -------------------------------------------------------------
# Pairs are matched without regard to order. exact agrees with A,
# distinct with C or E. Any other label counts against a definite claim,
# F included.

key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))
tp <- audit$pairs
lp$tool <- tp$relationship[match(key(lp$row_a, lp$row_b), key(tp$row_a, tp$row_b))]
lp$rule <- tp$rule_id[match(key(lp$row_a, lp$row_b), key(tp$row_a, tp$row_b))]
if (anyNA(lp$tool)) stop("Some labeled pairs are missing from the audit.")

definite <- lp$tool %in% c("exact", "distinct")
lp$c1 <- ifelse(!definite, NA,
         ifelse(lp$tool == "exact", lp$category == "A", lp$category %in% c("C", "E")))
c1_right <- sum(lp$c1, na.rm = TRUE)
c1_total <- sum(definite)


# ---- C2, groups ------------------------------------------------------------
# Group label: most severe pair label, E > C > F > D > B > A
# (notes/blind2_labels.md). Outcomes as in scripts/21_evaluate_v0.1.R.

lab_order <- c("E", "C", "F", "D", "B", "A")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")
gl <- tapply(lp$category, lp$group_id, function(x) lab_order[min(match(x, lab_order))])
groups <- data.frame(group_id = names(gl), category = unname(gl), stringsAsFactors = FALSE)
tg <- audit$groups
groups$tool <- tg$worst_relationship[match(groups$group_id, tg$group_id)]
groups$rule <- tg$worst_rule[match(groups$group_id, tg$group_id)]
if (anyNA(groups$tool)) stop("Some labeled groups are missing from the audit.")
groups$outcome <- with(groups, ifelse(category == "E", ifelse(tool == "exact", "E_exact", "E_not_exact"),
                               ifelse(tool == expected_map[category], "agree",
                               ifelse(tool == "unresolved", "safe_miss", "wrong"))))
in_c2 <- groups$category %in% c("A", "B", "C", "D")
c2_wrong <- sum(groups$outcome[in_c2] == "wrong")
c2_total <- sum(in_c2)


# ---- C3, lipid matches -----------------------------------------------------
# Every proposal across A and B that involves a sampled B row, set in the
# order (B row, A row). It is right when the labels hold that pair.

lm <- audit$lipid_matches
lm <- lm[lm$file_a != lm$file_b, ]
lm$b <- ifelse(lm$file_a == "B", lm$row_a, lm$row_b)
lm$a <- ifelse(lm$file_a == "B", lm$row_b, lm$row_a)
lm <- lm[lm$b %in% ll$row_b, ]

truth <- do.call(rbind, lapply(seq_len(nrow(ll)), function(i) {
  m <- trimws(strsplit(ll$matches_in_A[i], "[;,]")[[1]])
  m <- m[nzchar(m) & m != "none"]
  if (length(m) == 0) return(NULL)
  data.frame(b = ll$row_b[i], a = m, relationship = ll$relationship[i], stringsAsFactors = FALSE)
}))
lm$in_labels <- paste(lm$b, lm$a) %in% paste(truth$b, truth$a)
lm$label_relationship <- truth$relationship[match(paste(lm$b, lm$a), paste(truth$b, truth$a))]
truth$found <- paste(truth$b, truth$a) %in% paste(lm$b, lm$a)
c3_right <- sum(lm$in_labels)
c3_total <- nrow(lm)


# ---- Report ----------------------------------------------------------------

pct <- function(x, n) if (n == 0) NA_real_ else round(100 * x / n, 1)
summary <- data.frame(
  criterion = c("C1: definite claims that agree with the label (at least 90%)",
                "C2: groups labeled A to D with a wrong definite answer (at most 5%)",
                "C3: lipid proposals that are pairs in the labels (at least 90%)",
                "C3 recall: labeled lipid pairs the tool proposed (no threshold)"),
  count   = c(c1_right, c2_wrong, c3_right, sum(truth$found)),
  of      = c(c1_total, c2_total, c3_total, nrow(truth)),
  percent = c(pct(c1_right, c1_total), pct(c2_wrong, c2_total),
              pct(c3_right, c3_total), pct(sum(truth$found), nrow(truth))),
  stringsAsFactors = FALSE
)
summary$met <- c(if (c1_total == 0) NA else c1_right / c1_total >= 0.9,
                 c2_wrong / c2_total <= 0.05,
                 if (c3_total == 0) NA else c3_right / c3_total >= 0.9,
                 NA)

message("\n==== Second blind test")
print(summary, row.names = FALSE, right = FALSE)
message("\nPairs, label against tool:")
print(table(label = lp$category, tool = factor(lp$tool, levels = .severity)))
message("\nGroups, label against outcome:")
print(table(label = groups$category, outcome = groups$outcome))
message("\nLipid proposals involving sampled B rows, by tool relationship and label:")
print(table(tool = lm$relationship, label = ifelse(lm$in_labels, lm$label_relationship, "not in labels")))

# Counts only: committed.
write.csv(summary, file.path(pub, "summary.csv"), row.names = FALSE)
write.csv(as.data.frame(table(label = lp$category, tool = lp$tool, rule = lp$rule)) |>
            subset(Freq > 0), file.path(pub, "pairs_by_label_tool_rule.csv"), row.names = FALSE)
write.csv(as.data.frame(table(label = groups$category, outcome = groups$outcome)) |>
            subset(Freq > 0), file.path(pub, "groups_by_label_outcome.csv"), row.names = FALSE)
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
write.csv(truth, file.path(priv, "lipid_labels_found.csv"), row.names = FALSE, na = "")
write_audit(audit, file.path(priv, "audit"))
message("\nCounts in ", pub, "; details in ", priv, " (not committed).")
