# Acceptance checks for v0.7 (notes/decisions.md, 2026-10-07, v0.7 plan
# and the feasibility entry after it). All three parts are after
# adaptation.
#
#   Development data: evaluation against notes/review.csv, every pair
#     whose rule changed against v0.6.0, lipid match proposals lost and
#     gained.
#   COMETS: every pair whose rule changed, with the expert's view, and
#     definite claims against the expert.
#   Second blind test files: scored as in scripts/11_blind_test_2/04_test.R.
#
# Every run gets a copy of the PubChem cache of the v0.6.0 run it is
# compared with, so only the package differs.
# Counts go to results/v07_check/. Rows of the private file A go to
# data-raw/own/blind2/results_v07/.
#
# Run from the project root: Rscript scripts/12_v0.7/02_acceptance.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.7.0")
options(width = 220)

out  <- "results/v07_check"
priv <- "data-raw/own/blind2/results_v07"
dir.create(out,  showWarnings = FALSE, recursive = TRUE)
dir.create(priv, showWarnings = FALSE, recursive = TRUE)

lmsd <- readRDS("data-raw/lipidmaps/lmsd_db.rds")
stopifnot(identical(lmsd$md5, "134c2b5c6d76e9883d46339a6289eefd"))
sl <- load_swisslipids("data-raw/swisslipids/lipids.tsv", expected_md5 = "b6ddd5bc3003967dce95e905e5748953")

rerun <- function(old, cache_from, cache_to) {
  if (!file.exists(cache_to)) file.copy(cache_from, cache_to)
  audit_metlinkr(old$info$files, swisslipids = sl, cache_file = cache_to, lmsd = lmsd)
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

show <- function(ch, cols) {
  if (nrow(ch) == 0) return(invisible())
  print(ch[, cols], row.names = FALSE, right = FALSE)
}


# ---- Development data --------------------------------------------------------

old_dev <- readRDS("results/v06_check/audit_development.rds")
dev <- rerun(old_dev, "results/v06_check/pubchem_cache_development.rds",
             file.path(out, "pubchem_cache_development.rds"))
saveRDS(dev, file.path(out, "audit_development.rds"))

message("==== development")
message("Rows with an isomer marker: ", sum(dev$rows$name_isomer_marker))
ch <- changed_pairs(old_dev, dev)
message("Pairs whose rule changed against v0.6.0: ", nrow(ch))
show(ch, c("name_a", "name_b", "rule_id_old", "rule_id_new", "relationship_new"))
write.csv(ch, file.path(out, "dev_changed_pairs.csv"), row.names = FALSE)

review <- read.csv("notes/review.csv", colClasses = "character")
expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")
rel <- dev$groups$worst_relationship[match(review$metlinkr_name, dev$groups$group_id)]
outcome <- ifelse(review$category == "E", ifelse(rel == "exact", "E_exact", "E_not_exact"),
           ifelse(rel == unname(expected_map[review$category]), "agree",
           ifelse(rel == "unresolved", "safe_miss", "wrong")))
message("Evaluation against the review (v0.6.0: 35 agree, 8 safe_miss, 0 wrong, 10 E_not_exact):")
print(table(outcome))
ld <- lipid_diff(old_dev, dev)
message("Lipid match proposals: ", nrow(old_dev$lipid_matches), " before, ", nrow(dev$lipid_matches),
        " now; lost ", nrow(ld$lost), ", gained ", nrow(ld$gained))
show(ld$lost, c("name_a", "name_b", "relationship"))
show(ld$gained, c("name_a", "name_b", "relationship"))
write.csv(ld$lost, file.path(out, "dev_lipid_matches_lost.csv"), row.names = FALSE)
write.csv(ld$gained, file.path(out, "dev_lipid_matches_gained.csv"), row.names = FALSE)
dev_summary <- c(table(factor(outcome, levels = c("agree", "safe_miss", "wrong", "E_exact", "E_not_exact"))),
                 changed_pairs = nrow(ch), lipid_lost = nrow(ld$lost), lipid_gained = nrow(ld$gained))


# ---- COMETS ------------------------------------------------------------------

old_c <- readRDS("results/v06_check/audit_comets.rds")
com <- rerun(old_c, "results/v06_check/pubchem_cache_comets.rds", file.path(out, "pubchem_cache_comets.rds"))
saveRDS(com, file.path(out, "audit_comets.rds"))

# Expert groups, joined as in scripts/08_comets/02_blind_test.R and scripts/10_v0.6/04_acceptance.R.
src <- "data-raw/comets"
norm_id <- function(x) toupper(gsub("[[:space:]]", "", x))
files <- com$info$files
ids <- unlist(lapply(files$label, function(s)
  read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)[[1]]))
rk <- paste(com$rows$file, norm_id(ids), sep = "|")
mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character", na.strings = c("", "NA"))
mk <- paste(sub("\\.xlsx$", "", mm$SOURCE_FILE), norm_id(mm$IPT_METABID), sep = "|")
groups_of <- split(mm$Manual_Metabolite_Group_Label, mk)[rk]
names(groups_of) <- com$rows$row_id
same <- function(a, b) vapply(seq_along(a), function(i)
  length(intersect(groups_of[[a[i]]], groups_of[[b[i]]])) > 0, logical(1))

message("\n==== comets (after adaptation)")
message("Rows with an isomer marker: ", sum(com$rows$name_isomer_marker))
chc <- changed_pairs(old_c, com)
chc$expert_same <- if (nrow(chc) > 0) same(chc$row_a, chc$row_b) else logical()
message("Pairs whose rule changed against v0.6.0: ", nrow(chc))
show(chc, c("name_a", "name_b", "rule_id_old", "rule_id_new", "relationship_new", "expert_same"))
write.csv(chc, file.path(out, "comets_changed_pairs.csv"), row.names = FALSE)

# Definite claims on cross-file pairs, against the expert, for both runs.
definite_agree <- function(audit) {
  f <- audit$rows$file[match(audit$pairs$row_a, audit$rows$row_id)] !=
       audit$rows$file[match(audit$pairs$row_b, audit$rows$row_id)]
  d <- audit$pairs[f & audit$pairs$relationship %in% c("exact", "distinct"), ]
  s <- same(d$row_a, d$row_b)
  c(right = sum(d$relationship == "exact" & s) + sum(d$relationship == "distinct" & !s), total = nrow(d))
}
d06 <- definite_agree(old_c); d07 <- definite_agree(com)
right <- d07[["right"]]
message("Definite claims across files that agree with the expert: v0.6.0 ", d06[["right"]], " of ", d06[["total"]],
        ", v0.7.0 ", d07[["right"]], " of ", d07[["total"]])
lc <- lipid_diff(old_c, com)
message("Lipid match proposals: lost ", nrow(lc$lost), ", gained ", nrow(lc$gained))
comets_summary <- c(changed_pairs = nrow(chc), definite_right_v06 = d06[["right"]], definite_total_v06 = d06[["total"]],
                    definite_right = right, definite_total = d07[["total"]],
                    lipid_lost = nrow(lc$lost), lipid_gained = nrow(lc$gained))


# ---- Second blind test files -------------------------------------------------

old_b <- readRDS("data-raw/own/blind2/results/audit.rds")
bl <- rerun(old_b, "data-raw/own/blind2/results/pubchem_cache.rds", file.path(priv, "pubchem_cache.rds"))
saveRDS(bl, file.path(priv, "audit.rds"))

lab <- "data-raw/own/blind2"
md5 <- function(p) unname(tools::md5sum(p))
stopifnot(md5(file.path(lab, "label_pairs.csv"))  == "3f70bedbd9a0c5a1f21204571a37a75d",
          md5(file.path(lab, "label_lipids.csv")) == "bd170037bd8d528ddfc5d984fc5e9813")
lp0 <- read.csv(file.path(lab, "label_pairs.csv"),  colClasses = "character", na.strings = "")
ll  <- read.csv(file.path(lab, "label_lipids.csv"), colClasses = "character", na.strings = "")
lab_order <- c("E", "C", "F", "D", "B", "A")
key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))

score <- function(audit) {
  lp <- lp0
  i <- match(key(lp$row_a, lp$row_b), key(audit$pairs$row_a, audit$pairs$row_b))
  lp$tool <- audit$pairs$relationship[i]; lp$rule <- audit$pairs$rule_id[i]
  definite <- lp$tool %in% c("exact", "distinct")
  lp$c1 <- ifelse(!definite, NA, ifelse(lp$tool == "exact", lp$category == "A", lp$category %in% c("C", "E")))
  gl <- tapply(lp$category, lp$group_id, function(x) lab_order[min(match(x, lab_order))])
  g <- data.frame(group_id = names(gl), category = unname(gl), stringsAsFactors = FALSE)
  g$tool <- audit$groups$worst_relationship[match(g$group_id, audit$groups$group_id)]
  g$outcome <- with(g, ifelse(category == "E", ifelse(tool == "exact", "E_exact", "E_not_exact"),
                       ifelse(tool == expected_map[category], "agree",
                       ifelse(tool == "unresolved", "safe_miss", "wrong"))))
  in_c2 <- g$category %in% c("A", "B", "C", "D")
  lm <- audit$lipid_matches[audit$lipid_matches$file_a != audit$lipid_matches$file_b, ]
  lm$b <- ifelse(lm$file_a == "B", lm$row_a, lm$row_b)
  lm$a <- ifelse(lm$file_a == "B", lm$row_b, lm$row_a)
  lm <- lm[lm$b %in% ll$row_b, ]
  truth <- do.call(rbind, lapply(seq_len(nrow(ll)), function(k) {
    m <- trimws(strsplit(ll$matches_in_A[k], "[;,]")[[1]]); m <- m[nzchar(m) & m != "none"]
    if (length(m) == 0) NULL else data.frame(b = ll$row_b[k], a = m, stringsAsFactors = FALSE)
  }))
  list(pairs = lp, groups = g,
       counts = c(c1_right = sum(lp$c1, na.rm = TRUE), c1_total = sum(definite),
                  c2_wrong = sum(g$outcome[in_c2] == "wrong"), c2_total = sum(in_c2),
                  c3_right = sum(paste(lm$b, lm$a) %in% paste(truth$b, truth$a)), c3_total = nrow(lm),
                  recall = sum(paste(truth$b, truth$a) %in% paste(lm$b, lm$a)), recall_total = nrow(truth)))
}
s06 <- score(old_b); s07 <- score(bl)
stopifnot(identical(unname(s06$counts[c("c1_right", "c1_total", "c2_wrong", "c2_total")]), c(48L, 65L, 12L, 58L)))

message("\n==== second blind test files (after adaptation, no threshold)")
message("Rows with an isomer marker: ", sum(bl$rows$name_isomer_marker), " (A: ",
        sum(bl$rows$name_isomer_marker & bl$rows$file == "A"), ")")
bcounts <- data.frame(measure = names(s07$counts), v0.6.0 = unname(s06$counts), v0.7.0 = unname(s07$counts))
print(bcounts, row.names = FALSE, right = FALSE)
pch <- merge(s06$pairs[, c("row_a", "row_b", "group_id", "category", "tool", "rule")],
             s07$pairs[, c("row_a", "row_b", "tool", "rule")], by = c("row_a", "row_b"), suffixes = c("_v06", "_v07"))
pch <- pch[pch$rule_v06 != pch$rule_v07, ]
message("Labeled pairs whose rule changed, by label:")
if (nrow(pch) > 0) print(as.data.frame(table(label = pch$category, v06 = paste(pch$tool_v06, pch$rule_v06),
                                              v07 = paste(pch$tool_v07, pch$rule_v07))) |> subset(Freq > 0),
                         row.names = FALSE, right = FALSE)
message("Groups, label against outcome:"); print(table(label = s07$groups$category, outcome = s07$groups$outcome))
ball <- changed_pairs(old_b, bl)
message("All pairs in the two files whose rule changed: ", nrow(ball), "; by change:")
print(table(paste(ball$rule_id_old, "->", ball$rule_id_new)))
write.csv(pch,  file.path(priv, "labeled_changed_pairs.csv"), row.names = FALSE, na = "")
write.csv(ball, file.path(priv, "all_changed_pairs.csv"), row.names = FALSE, na = "")


# ---- Counts ------------------------------------------------------------------

write.csv(data.frame(measure = names(dev_summary), value = unname(dev_summary)),
          file.path(out, "development_summary.csv"), row.names = FALSE)
write.csv(data.frame(measure = names(comets_summary), value = unname(comets_summary)),
          file.path(out, "comets_summary.csv"), row.names = FALSE)
write.csv(bcounts, file.path(out, "blind2_counts_after_adaptation.csv"), row.names = FALSE)
write.csv(as.data.frame(table(change = paste(ball$rule_id_old, "->", ball$rule_id_new))),
          file.path(out, "blind2_changed_rules.csv"), row.names = FALSE)
message("\nCounts in ", out, "; rows of file A in ", priv, ".")
