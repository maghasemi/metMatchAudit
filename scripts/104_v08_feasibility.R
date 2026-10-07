# v0.8 feasibility (notes/decisions.md, 2026-10-07, v0.8 plan, step 2).
# Changes nothing in the package. On the v0.7.0 audits of the four sets,
# lists every R09 and R09b pair where neither name states stereo, with
# the label it carries, and checks the go-ahead in the plan: at least two
# thirds of the labeled touched pairs have a label other than D.
#
# Rows of the private file A are listed only under
# data-raw/own/v08_feasibility/; counts go to results/v08_feasibility/.
#
# Run from the project root: Rscript scripts/104_v08_feasibility.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.7.0")
options(width = 220)

pub  <- "results/v08_feasibility"
priv <- "data-raw/own/v08_feasibility"
dir.create(pub,  showWarnings = FALSE, recursive = TRUE)
dir.create(priv, showWarnings = FALSE, recursive = TRUE)

touched <- function(audit) {
  p <- audit$pairs[audit$pairs$rule_id %in% c("R09", "R09b"), ]
  st <- stats::setNames(audit$rows$name_states_stereo, audit$rows$row_id)
  nm <- stats::setNames(audit$rows$original_name, audit$rows$row_id)
  p <- p[!st[p$row_a] & !st[p$row_b], ]
  p$name_a <- unname(nm[p$row_a]); p$name_b <- unname(nm[p$row_b])
  p
}
key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))
all_r09 <- function(audit) sum(audit$pairs$rule_id %in% c("R09", "R09b"))

# Development: label is the review category of the group.
dev <- readRDS("results/v07_check/audit_development.rds")
td <- touched(dev)
review <- read.csv("notes/review.csv", colClasses = "character")
td$label <- review$category[match(td$group_id, review$metlinkr_name)]
td$set <- "development"

# COMETS: whether the expert groups the pair (as in scripts/96_v07_check.R).
com <- readRDS("results/v07_check/audit_comets.rds")
tc <- touched(com)
src <- "data-raw/comets"
norm_id <- function(x) toupper(gsub("[[:space:]]", "", x))
ids <- unlist(lapply(com$info$files$label, function(s)
  read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)[[1]]))
rk <- paste(com$rows$file, norm_id(ids), sep = "|")
mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character", na.strings = c("", "NA"))
mk <- paste(sub("\\.xlsx$", "", mm$SOURCE_FILE), norm_id(mm$IPT_METABID), sep = "|")
groups_of <- split(mm$Manual_Metabolite_Group_Label, mk)[rk]
names(groups_of) <- com$rows$row_id
tc$expert_same <- vapply(seq_len(nrow(tc)), function(i)
  length(intersect(groups_of[[tc$row_a[i]]], groups_of[[tc$row_b[i]]])) > 0, logical(1))
tc$label <- NA_character_; tc$set <- "comets"

# Blind tests: pair labels.
pair_label <- function(t, sheet) {
  lp <- read.csv(sheet, colClasses = "character", na.strings = "")
  lp$category[match(key(t$row_a, t$row_b), key(lp$row_a, lp$row_b))]
}
b2 <- readRDS("data-raw/own/blind2/results_v07/audit.rds")
t2 <- touched(b2); t2$label <- pair_label(t2, "data-raw/own/blind2/label_pairs.csv"); t2$set <- "blind2"
b3 <- readRDS("data-raw/own/blind3/results/audit.rds")
t3 <- touched(b3); t3$label <- pair_label(t3, "data-raw/own/blind3/label_pairs.csv"); t3$set <- "blind3"

cols <- c("set", "group_id", "row_a", "row_b", "name_a", "name_b", "rule_id", "label")
counts <- data.frame(
  set = c("development", "comets", "blind2", "blind3"),
  r09_pairs = c(all_r09(dev), all_r09(com), all_r09(b2), all_r09(b3)),
  touched = c(nrow(td), nrow(tc), nrow(t2), nrow(t3)),
  labeled = c(sum(!is.na(td$label)), NA, sum(!is.na(t2$label)), sum(!is.na(t3$label))),
  label_not_D = c(sum(td$label != "D", na.rm = TRUE), NA, sum(t2$label != "D", na.rm = TRUE),
                  sum(t3$label != "D", na.rm = TRUE)),
  expert_same = c(NA, sum(tc$expert_same), NA, NA))
message("==== R09 and R09b pairs where neither name states stereo")
print(counts, row.names = FALSE)

lab <- rbind(td[, cols], t2[, cols], t3[, cols])
lab <- lab[!is.na(lab$label), ]
message("\nLabels of touched pairs (development, blind2, blind3):")
print(table(set = lab$set, label = lab$label))
share <- mean(lab$label != "D")
message(sprintf("\nShare with a label other than D: %d of %d (%.1f%%). Go-ahead needs at least 66.7%%: %s",
                sum(lab$label != "D"), nrow(lab), 100 * share, if (share >= 2 / 3) "MET" else "NOT MET"))

message("\nDevelopment and COMETS pairs (public data):")
print(rbind(td[, cols], transform(tc[, setdiff(cols, "label")], label = ifelse(tc$expert_same, "expert_same", "expert_diff"))),
      row.names = FALSE, right = FALSE)

write.csv(counts, file.path(pub, "counts.csv"), row.names = FALSE)
write.csv(as.data.frame(table(set = lab$set, label = lab$label)), file.path(pub, "labels_of_touched.csv"), row.names = FALSE)
write.csv(rbind(td[, cols], transform(tc[, setdiff(cols, "label")], label = ifelse(tc$expert_same, "expert_same", "expert_diff"))),
          file.path(pub, "touched_development_comets.csv"), row.names = FALSE)
write.csv(rbind(t2[, cols], t3[, cols]), file.path(priv, "touched_blind.csv"), row.names = FALSE)
message("\nCounts in ", pub, "; blind-test rows in ", priv, ".")
