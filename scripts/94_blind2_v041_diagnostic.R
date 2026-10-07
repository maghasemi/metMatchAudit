# Diagnostic after the second blind test: v0.4.1, the version frozen for
# the COMETS test, on the same two files and labels as v0.6.0 in
# scripts/93_blind2_test.R.
#
# This is not a blind result. It answers one question: did the changes
# after v0.4.1 make things worse on data neither version was built on?
# Both versions are scored with the same code below, and the v0.6.0 run
# must reproduce the numbers 93 reported.
#
# v0.4.1 is loaded from a temporary git worktree at its tag and gets a
# copy of the PubChem cache from the v0.6.0 run, so both versions see the
# same PubChem answers. v0.4.1 has no LMSD and no vendor table.
#
# Privacy as in 93: counts go to results/blind2_v041_diagnostic/, rows
# to data-raw/own/blind2/results_v041/.
#
# Run from the project root: Rscript scripts/94_blind2_v041_diagnostic.R

options(width = 200)
own  <- "data-raw/own"
lab  <- file.path(own, "blind2")
priv <- file.path(lab, "results_v041")
pub  <- "results/blind2_v041_diagnostic"
dir.create(priv, showWarnings = FALSE, recursive = TRUE)
dir.create(pub,  showWarnings = FALSE, recursive = TRUE)

md5 <- function(p) unname(tools::md5sum(p))
stopifnot(md5(file.path(lab, "label_pairs.csv"))  == "3f70bedbd9a0c5a1f21204571a37a75d",
          md5(file.path(lab, "label_lipids.csv")) == "bd170037bd8d528ddfc5d984fc5e9813")
lp0 <- read.csv(file.path(lab, "label_pairs.csv"),  colClasses = "character", na.strings = "")
ll  <- read.csv(file.path(lab, "label_lipids.csv"), colClasses = "character", na.strings = "")


# ---- v0.4.1 run ------------------------------------------------------------

v06 <- readRDS(file.path(lab, "results", "audit.rds"))
stopifnot(v06$info$package_version == "0.6.0")

v041_rds <- file.path(priv, "audit.rds")
if (!file.exists(v041_rds)) {
  wt <- file.path(tempdir(), "mma_v0.4.1")
  if (system(paste("git worktree add --detach", shQuote(wt), "v0.4.1")) != 0) stop("git worktree failed.")
  tryCatch({
    pkgload::load_all(wt, quiet = TRUE, export_all = FALSE)
    stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.4.1")

    cache <- file.path(priv, "pubchem_cache.rds")
    if (!file.exists(cache)) file.copy(file.path(lab, "results", "pubchem_cache.rds"), cache)

    # The files table of the v0.6.0 run, which has no vendor column.
    a41 <- audit_metlinkr(v06$info$files,
                          swisslipids  = "data-raw/swisslipids/lipids.tsv",
                          expected_md5 = "b6ddd5bc3003967dce95e905e5748953",
                          cache_file   = cache)
    saveRDS(a41, v041_rds)
  }, finally = system(paste("git worktree remove --force", shQuote(wt))))
}
a41 <- readRDS(v041_rds)
stopifnot(a41$info$package_version == "0.4.1")


# ---- Scoring, as in scripts/93_blind2_test.R -------------------------------

severity  <- c("distinct", "unresolved", "broader_narrower", "stereo_variant", "exact")
lab_order <- c("E", "C", "F", "D", "B", "A")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")
key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))

score <- function(audit) {
  lp <- lp0
  tp <- audit$pairs
  i <- match(key(lp$row_a, lp$row_b), key(tp$row_a, tp$row_b))
  if (anyNA(i)) stop("Labeled pairs missing from the audit.")
  lp$tool <- tp$relationship[i]; lp$rule <- tp$rule_id[i]
  definite <- lp$tool %in% c("exact", "distinct")
  lp$c1 <- ifelse(!definite, NA,
           ifelse(lp$tool == "exact", lp$category == "A", lp$category %in% c("C", "E")))

  gl <- tapply(lp$category, lp$group_id, function(x) lab_order[min(match(x, lab_order))])
  g <- data.frame(group_id = names(gl), category = unname(gl), stringsAsFactors = FALSE)
  j <- match(g$group_id, audit$groups$group_id)
  g$tool <- audit$groups$worst_relationship[j]; g$rule <- audit$groups$worst_rule[j]
  g$outcome <- with(g, ifelse(category == "E", ifelse(tool == "exact", "E_exact", "E_not_exact"),
                       ifelse(tool == expected_map[category], "agree",
                       ifelse(tool == "unresolved", "safe_miss", "wrong"))))
  in_c2 <- g$category %in% c("A", "B", "C", "D")

  lm <- audit$lipid_matches
  lm <- lm[lm$file_a != lm$file_b, ]
  lm$b <- ifelse(lm$file_a == "B", lm$row_a, lm$row_b)
  lm$a <- ifelse(lm$file_a == "B", lm$row_b, lm$row_a)
  lm <- lm[lm$b %in% ll$row_b, ]
  truth <- do.call(rbind, lapply(seq_len(nrow(ll)), function(k) {
    m <- trimws(strsplit(ll$matches_in_A[k], "[;,]")[[1]])
    m <- m[nzchar(m) & m != "none"]
    if (length(m) == 0) return(NULL)
    data.frame(b = ll$row_b[k], a = m, stringsAsFactors = FALSE)
  }))
  lm$in_labels <- paste(lm$b, lm$a) %in% paste(truth$b, truth$a)

  list(pairs = lp, groups = g,
       counts = c(c1_right = sum(lp$c1, na.rm = TRUE), c1_total = sum(definite),
                  c2_wrong = sum(g$outcome[in_c2] == "wrong"), c2_total = sum(in_c2),
                  c3_right = sum(lm$in_labels), c3_total = nrow(lm),
                  recall = sum(paste(truth$b, truth$a) %in% paste(lm$b, lm$a)), recall_total = nrow(truth),
                  lipid_matches_all = nrow(audit$lipid_matches)))
}

s06 <- score(v06)
s41 <- score(a41)

# The v0.6.0 numbers must match what 93 reported.
stopifnot(identical(unname(s06$counts[c("c1_right", "c1_total", "c2_wrong", "c2_total")]),
                    c(48L, 65L, 12L, 58L)))


# ---- Report ----------------------------------------------------------------

counts <- data.frame(measure = names(s41$counts), v0.4.1 = unname(s41$counts),
                     v0.6.0 = unname(s06$counts), stringsAsFactors = FALSE)
message("\n==== v0.4.1 against v0.6.0 on the second blind test files (diagnostic, not blind)")
print(counts, row.names = FALSE, right = FALSE)

ch <- merge(s41$pairs[, c("group_id", "row_a", "row_b", "category", "tool", "rule", "c1")],
            s06$pairs[, c("row_a", "row_b", "tool", "rule", "c1")],
            by = c("row_a", "row_b"), suffixes = c("_v041", "_v06"))
ch <- ch[ch$rule_v041 != ch$rule_v06, ]
message("\nLabeled pairs whose rule differs: ", nrow(ch))
if (nrow(ch) > 0) print(as.data.frame(table(label = ch$category,
                                            v041 = paste(ch$tool_v041, ch$rule_v041),
                                            v06 = paste(ch$tool_v06, ch$rule_v06))) |> subset(Freq > 0),
                        row.names = FALSE, right = FALSE)

gch <- merge(s41$groups[, c("group_id", "category", "outcome")], s06$groups[, c("group_id", "outcome")],
             by = "group_id", suffixes = c("_v041", "_v06"))
gch <- gch[gch$outcome_v041 != gch$outcome_v06, ]
message("\nGroups whose outcome differs: ", nrow(gch))
if (nrow(gch) > 0) print(table(label = gch$category, change = paste(gch$outcome_v041, "->", gch$outcome_v06)))

pairs_tab <- as.data.frame(table(label = ch$category, v041 = paste(ch$tool_v041, ch$rule_v041),
                                 v06 = paste(ch$tool_v06, ch$rule_v06)))
write.csv(counts, file.path(pub, "counts.csv"), row.names = FALSE)
write.csv(subset(pairs_tab, Freq > 0), file.path(pub, "changed_pairs_by_label.csv"), row.names = FALSE)
write.csv(as.data.frame(table(label = gch$category, v041 = gch$outcome_v041, v06 = gch$outcome_v06)) |>
            subset(Freq > 0), file.path(pub, "changed_groups_by_label.csv"), row.names = FALSE)
write.csv(ch,  file.path(priv, "changed_pairs.csv"),  row.names = FALSE, na = "")
write.csv(gch, file.path(priv, "changed_groups.csv"), row.names = FALSE, na = "")
message("\nCounts in ", pub, "; details in ", priv, " (not committed).")
