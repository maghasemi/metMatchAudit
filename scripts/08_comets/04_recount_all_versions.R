# Recounts the COMETS comparison for every version on disk, with the code
# of scripts/08_comets/02_blind_test.R. The README gives 123 of 123 definite
# claims for v0.6, while scripts/12_v0.7/02_acceptance.R counts 122 of 122 for
# both v0.6.0 and v0.7.0. This script finds out which number is right.
#
# The v0.4.1 run must reproduce the blind result (121 of 125, F1 0.893
# against 0.895); otherwise the script stops, since the recount would
# not be the same computation.
#
# Output: results/comets_recount/. Run from the project root:
# Rscript scripts/08_comets/04_recount_all_versions.R

options(width = 200)
src <- "data-raw/comets"
out <- "results/comets_recount"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

runs <- c(v0.4.1 = "results/comets_test/audit.rds",
          v0.4.2 = "results/comets_v0.4.2/audit.rds",
          v0.6.0 = "results/v06_check/audit_comets.rds",
          v0.7.0 = "results/v07_check/audit_comets.rds")

norm_id <- function(x) toupper(gsub("[[:space:]]", "", x))
mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character", na.strings = c("", "NA"))
mm_key <- paste(sub("\\.xlsx$", "", mm$SOURCE_FILE), norm_id(mm$IPT_METABID), sep = "|")

prf <- function(tp, predicted, positives) {
  p <- if (predicted == 0) NA else tp / predicted
  r <- tp / positives
  f <- if (is.na(p) || (p + r) == 0) NA else 2 * p * r / (p + r)
  c(precision = round(p, 3), recall = round(r, 3), f1 = round(f, 3))
}

evaluate <- function(audit) {
  rows <- audit$rows
  labels <- audit$info$files$label
  ids <- unlist(lapply(labels, function(s)
    read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)[[1]]))
  stopifnot(length(ids) == nrow(rows))
  row_key <- paste(rows$file, norm_id(ids), sep = "|")
  membership <- unique(data.frame(row_id = rows$row_id[match(mm_key, row_key)],
                                  group = mm$Manual_Metabolite_Group_Label, stringsAsFactors = FALSE))
  stopifnot(!anyNA(membership$row_id))
  groups_of <- split(membership$group, membership$row_id)
  file_of <- stats::setNames(rows$file, rows$row_id)
  same <- function(a, b) vapply(seq_along(a), function(i)
    length(intersect(groups_of[[a[i]]], groups_of[[b[i]]])) > 0, logical(1))

  per_group <- split(membership$row_id, membership$group)
  pp <- do.call(cbind, lapply(per_group, function(r) {
    if (length(r) < 2) return(NULL)
    cm <- utils::combn(sort(r), 2)
    cm[, file_of[cm[1, ]] != file_of[cm[2, ]], drop = FALSE]
  }))
  positives <- length(unique(paste(pp[1, ], pp[2, ], sep = "\r")))

  p <- audit$pairs
  p <- p[file_of[p$row_a] != file_of[p$row_b], ]
  p$expert_same <- same(p$row_a, p$row_b)
  ex <- p[p$relationship == "exact", ]; di <- p[p$relationship == "distinct", ]
  alone <- prf(sum(p$expert_same), nrow(p), positives)
  kept <- p[p$relationship != "distinct", ]
  with  <- prf(sum(kept$expert_same), nrow(kept), positives)
  list(pairs = p,
       counts = c(cross_pairs = nrow(p), exact = nrow(ex), exact_right = sum(ex$expert_same),
                  distinct = nrow(di), distinct_right = sum(!di$expert_same),
                  definite_right = sum(ex$expert_same) + sum(!di$expert_same),
                  definite_total = nrow(ex) + nrow(di),
                  precision_alone = alone[["precision"]], f1_alone = alone[["f1"]],
                  precision_with = with[["precision"]], f1_with = with[["f1"]]))
}

res <- lapply(runs, function(f) evaluate(readRDS(f)))
tab <- do.call(cbind, lapply(res, `[[`, "counts"))
message("==== COMETS, recounted with the code of scripts/08_comets/02_blind_test.R")
print(tab)

stopifnot(tab["definite_right", "v0.4.1"] == 121, tab["definite_total", "v0.4.1"] == 125,
          tab["f1_with", "v0.4.1"] == 0.893, tab["f1_alone", "v0.4.1"] == 0.895)

# Pairs whose relationship differs from one version to the next.
key <- function(p) paste(p$row_a, p$row_b)
for (i in 2:length(res)) {
  a <- res[[i - 1]]$pairs; b <- res[[i]]$pairs
  m <- merge(a[, c("row_a", "row_b", "relationship", "rule_id", "expert_same")],
             b[, c("row_a", "row_b", "relationship", "rule_id")],
             by = c("row_a", "row_b"), suffixes = c("_old", "_new"))
  m <- m[m$relationship_old != m$relationship_new, ]
  message("\n", names(res)[i - 1], " -> ", names(res)[i], ": ", nrow(m), " cross-file pairs changed relationship")
  if (nrow(m) > 0) print(as.data.frame(table(old = m$relationship_old, new = m$relationship_new,
                                             expert_same = m$expert_same)) |> subset(Freq > 0),
                         row.names = FALSE)
}
write.csv(data.frame(measure = rownames(tab), tab, check.names = FALSE),
          file.path(out, "counts_by_version.csv"), row.names = FALSE)
