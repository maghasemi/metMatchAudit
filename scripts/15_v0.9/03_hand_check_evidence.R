# Evidence for the v0.9 hand check, and the sheet filled from it
# (notes/decisions.md, 2026-11-07, who did the v0.9 hand check).
#
# For each of the 23 pairs in results/v09_check/hand_check.csv:
#   1. every alternative in the Metabolon name must add up to the chain sum
#      in the Broad name (a plasmalogen "P-" chain counts one double bond);
#   2. the formula of the Metabolon SMILES (RDKit) must equal the formula
#      worked out from the Broad name alone, from the lipid class;
#   3. for a claim, the Metabolon row must be the one with more detail.
#
# Neither check uses the parser or the rules of the package. Broad gives no
# m/z and no formula for these rows, so this compares names and structures,
# not measurements.
#
# Run from the project root: Rscript scripts/15_v0.9/03_hand_check_evidence.R

pkgload::load_all(".", quiet = TRUE)
options(width = 200)

sheet <- "results/v09_check/hand_check.csv"
h <- read.csv(sheet, colClasses = "character", na.strings = NULL)
audit <- readRDS("results/v09_check/audit_development.rds")
files <- audit$info$files
metabolon <- data.table::fread(files$path[files$label == "LEOCC_Metabolon"], colClasses = "character",
                               encoding = "Latin-1", na.strings = c("", "NA"), data.table = FALSE)

u <- h[!duplicated(h[, c("row_a", "row_b")]), c("row_a", "row_b", "name_a", "name_b", "rule", "relationship")]
u$smiles <- metabolon$SMILES[as.integer(sub(".*:", "", u$row_a))]

# Neutral formula from the Broad name, by lipid class.
class_formula <- function(name) {
  g <- regmatches(name, regexec("^(LPC O-|LPC |PC O-|PC |PE O-|SM )([0-9]+):([0-9]+)", name))[[1]]
  n <- as.integer(g[3]); d <- as.integer(g[4])
  p <- switch(trimws(g[2]),
    "LPC"    = c(n + 8, 2 * n + 18 - 2 * d, 1, 7),
    "LPC O-" = c(n + 8, 2 * n + 20 - 2 * d, 1, 6),
    "PC"     = c(n + 8, 2 * n + 16 - 2 * d, 1, 8),
    "PC O-"  = c(n + 8, 2 * n + 18 - 2 * d, 1, 7),
    "PE O-"  = c(n + 5, 2 * n + 12 - 2 * d, 1, 7),
    "SM"     = c(n + 5, 2 * n + 13 - 2 * d, 2, 6))
  sprintf("C%dH%d%sO%dP", p[1], p[2], if (p[3] == 1) "N" else paste0("N", p[3]), p[4])
}

# Chain sum of each alternative in the last brackets of a Metabolon name.
chain_sums <- function(name) {
  inner <- sub(".*\\(([^()]*)\\)[^()]*$", "\\1", name)
  vapply(trimws(strsplit(inner, ",")[[1]]), function(alt) {
    chains <- regmatches(alt, gregexpr("[dP]?-?[0-9]+:[0-9]+", alt))[[1]]
    v <- do.call(rbind, strsplit(sub("^[dP]?-?", "", chains), ":"))
    paste0(sum(as.integer(v[, 1])), ":", sum(as.integer(v[, 2])) + sum(grepl("^P-", chains)))
  }, "", USE.NAMES = FALSE)
}

mods <- metMatchAudit:::.rdkit()
desc <- reticulate::import("rdkit.Chem.rdMolDescriptors")
u$formula_smiles <- vapply(u$smiles, function(s) desc$CalcMolFormula(mods$Chem$MolFromSmiles(s)), "", USE.NAMES = FALSE)
u$formula_name <- vapply(u$name_b, class_formula, "", USE.NAMES = FALSE)
u$formula_ok <- u$formula_smiles == u$formula_name

sums <- lapply(u$name_a, chain_sums)
u$alternatives <- lengths(sums)
u$sum_b <- sub("^[^0-9]*([0-9]+:[0-9]+).*$", "\\1", u$name_b)
u$sum_ok <- mapply(function(s, b) all(s == b), sums, u$sum_b)

level <- stats::setNames(audit$rows$lipid_level, audit$rows$row_id)
u$level_a <- unname(level[u$row_a])
u$level_b <- unname(level[u$row_b])
claim <- u$relationship != "unresolved"
u$detail_ok <- !claim | (u$level_a %in% c("structure", "group") & u$level_b %in% c("sn_position", "species"))

u$right <- ifelse(u$sum_ok & u$formula_ok & u$detail_ok, "yes", "no")
u$note <- ifelse(u$rule == "R99", "Same lipid by chain sum and formula. Safe miss: both rows are at structure level, and no rule covers that without two keys.",
          ifelse(u$rule == "R08", "All alternatives add up to the Species and the formula matches. Safe miss: the group row has no lipid profile.",
          ifelse(u$level_a == "structure" & u$level_b == "sn_position", "Broad gives no double bond position, so broader_narrower can be defended. A curator might call this the same compound.", "")))

print(u[, c("name_a", "name_b", "rule", "alternatives", "sum_ok", "formula_smiles", "formula_ok", "level_a", "level_b", "detail_ok", "right")],
      row.names = FALSE, right = FALSE)
message("Pairs: ", nrow(u), "; chain sums right: ", sum(u$sum_ok), "; formulas right: ", sum(u$formula_ok))
message("Claims: ", sum(claim), "; right: ", sum(u$right[claim] == "yes"), ". Unresolved, reported apart: ", sum(!claim))

write.csv(u, "results/v09_check/hand_check_evidence.csv", row.names = FALSE)
i <- match(paste(h$row_a, h$row_b), paste(u$row_a, u$row_b))
h$right <- u$right[i]
h$note <- u$note[i]
write.csv(h, sheet, row.names = FALSE, na = "")
