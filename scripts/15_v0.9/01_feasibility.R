# v0.9 feasibility (notes/decisions.md, 2026-11-07, v0.9 plan, step 2).
# Changes nothing in the package. Checks whether column `name` of the
# Broad example file can stand in for an empty `Metabolite` cell.
#
# Counts go to results/v09_feasibility/. Rows of the private file A are
# only counted, never listed.
#
# Run from the project root: Rscript scripts/15_v0.9/01_feasibility.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.7.1")
options(width = 200)

out <- "results/v09_feasibility"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

dev <- readRDS("results/v07_check/audit_development.rds")
path <- dev$info$files$path[dev$info$files$label == "Broad_2022Aug"]
d <- data.table::fread(path, colClasses = "character", encoding = "Latin-1",
                       na.strings = c("", "NA"), data.table = FALSE, showProgress = FALSE)
d$row_id <- paste0("Broad_2022Aug:", seq_len(nrow(d)))
empty <- function(x) is.na(x) | !nzchar(trimws(x))
blank <- empty(d$Metabolite)
message("Broad rows: ", nrow(d), "; empty Metabolite: ", sum(blank),
        "; of those with a value in `name`: ", sum(blank & !empty(d$name)))

# 1. Do the names of the blank rows parse?
p <- parse_lipid_names(d$name[blank])
share <- mean(p$parsed)
message("\n1. Blank rows whose `name` parses: ", sum(p$parsed), " of ", nrow(p),
        " (", round(100 * share, 1), "%). Needed: at least 80%.")
print(table(level = p$level, useNA = "ifany"))
if (any(!p$parsed)) {
  message("Names that did not parse:")
  print(d$name[blank][!p$parsed])
}

# 2. Where both columns are filled, do they mean the same lipid?
both <- !blank & !empty(d$name)
message("\n2. Rows with both columns filled: ", sum(both),
        "; identical text: ", sum(d$Metabolite[both] == d$name[both]))
pm <- parse_lipid_names(d$Metabolite[both])
pn <- parse_lipid_names(d$name[both])
ok <- pm$parsed & pn$parsed
same <- pm$species[ok] == pn$species[ok]
agree <- if (any(ok)) mean(same) else NA_real_
message("Both parse: ", sum(ok), "; same Species: ", sum(same),
        " (", round(100 * agree, 1), "%). Needed: at least 95%.")
if (any(!same)) {
  print(data.frame(Metabolite = d$Metabolite[both][ok][!same], name = d$name[both][ok][!same],
                   species_metabolite = pm$species[ok][!same], species_name = pn$species[ok][!same]),
        row.names = FALSE, right = FALSE)
}

# 3. Blank names in every file of the four test sets.
sets <- c(development = "results/v07_check/audit_development.rds",
          comets      = "results/v07_check/audit_comets.rds",
          blind2      = "data-raw/own/blind2/results_v07/audit.rds",
          blind3      = "data-raw/own/blind3/results/audit.rds")
counts <- do.call(rbind, lapply(names(sets), function(s) {
  if (!file.exists(sets[[s]])) {
    message("Not found, skipped: ", sets[[s]])
    return(NULL)
  }
  r <- readRDS(sets[[s]])$rows
  tab <- tapply(empty(r$original_name), r$file, sum)
  data.frame(set = s, file = names(tab), rows = as.vector(table(r$file)[names(tab)]),
             blank_names = as.vector(tab))
}))
message("\n3. Blank names per file:")
print(counts, row.names = FALSE, right = FALSE)
write.csv(counts, file.path(out, "blank_names.csv"), row.names = FALSE)

# 4. Development pairs that hold a blank row.
ids <- d$row_id[blank]
pr <- dev$pairs[dev$pairs$row_a %in% ids | dev$pairs$row_b %in% ids, ]
nm <- stats::setNames(dev$rows$original_name, dev$rows$row_id)
fill <- stats::setNames(d$name, d$row_id)
label <- function(id) ifelse(id %in% ids, paste0(fill[id], " [from name]"), nm[id])
pr$name_a <- unname(label(pr$row_a))
pr$name_b <- unname(label(pr$row_b))
message("\n4. Development pairs with a blank row: ", nrow(pr))
print(pr[, c("name_a", "name_b", "relationship", "rule_id")], row.names = FALSE, right = FALSE)
write.csv(pr[, c("row_a", "row_b", "group_id", "name_a", "name_b", "relationship", "rule_id")],
          file.path(out, "dev_pairs_with_blank_row.csv"), row.names = FALSE)

go <- share >= 0.8 && isTRUE(agree >= 0.95)
write.csv(data.frame(measure = c("blank_rows", "parsed", "parsed_share", "both_parse", "same_species", "same_share", "go_ahead"),
                     value = c(sum(blank), sum(p$parsed), round(share, 4), sum(ok), sum(same), round(agree, 4), go)),
          file.path(out, "summary.csv"), row.names = FALSE)
message("\nGo-ahead: ", if (go) "YES" else "NO")