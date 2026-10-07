# Second blind test, step 4: draws the groups and lipid rows to label and
# writes the two labeling sheets. Seed and rules are in notes/decisions.md
# (2026-10-07, second blind test). Reads only MetLinkR's output; the audit
# is not run here.
#
# The sheets hold names from the private file A, so they stay in
# data-raw/own/blind2/. Their md5 goes into notes/blind2_labels.md once
# they are filled in, which fixes the labels in git without publishing
# them.
#
# Run from the project root: Rscript scripts/11_blind_test_2/03_sample.R

pkgload::load_all(".", quiet = TRUE)

dir  <- "data-raw/own"
mout <- file.path(dir, "metlinkr_output")
out  <- file.path(dir, "blind2")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (file.exists(file.path(out, "label_pairs.csv"))) stop("Sheets exist already. Delete them on purpose to redraw.")

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
rows <- read_metlinkr(files)
message("Rows: ", sum(rows$file == "A"), " in A, ", sum(rows$file == "B"), " in B")

# Groups with rows from both files, in a fixed order before sampling.
g <- split(rows$file, rows$group_id)
both <- sort(names(g)[vapply(g, function(f) all(c("A", "B") %in% f), logical(1))])
message("MetLinkR groups with rows from both files: ", length(both))

set.seed(20261007)
picked <- if (length(both) <= 60) both else sample(both, 60)

pairs <- make_pairs(rows[rows$group_id %in% picked, ])
r <- function(id, col) rows[[col]][match(id, rows$row_id)]
sheet <- data.frame(
  group_id = pairs$group_id,
  row_a = pairs$row_a, name_a = r(pairs$row_a, "original_name"),
  pubchem_a = r(pairs$row_a, "pubchem"), hmdb_a = r(pairs$row_a, "hmdb"),
  inchikey_a = r(pairs$row_a, "inchikey"),
  row_b = pairs$row_b, name_b = r(pairs$row_b, "original_name"),
  pubchem_b = r(pairs$row_b, "pubchem"), hmdb_b = r(pairs$row_b, "hmdb"),
  inchikey_b = r(pairs$row_b, "inchikey"),
  cross_file = r(pairs$row_a, "file") != r(pairs$row_b, "file"),
  category = "", evidence = "",
  stringsAsFactors = FALSE
)
sheet <- sheet[order(sheet$group_id, !sheet$cross_file), ]
write.csv(sheet, file.path(out, "label_pairs.csv"), row.names = FALSE, na = "")
message("Groups drawn: ", length(picked), "; pairs to label: ", nrow(sheet),
        " (", sum(sheet$cross_file), " across files)")

# Lipid rows of B, by the class pattern fixed in the plan.
lip_pat <- "^(L?P[ACEGIS]|SM|Cer|HexCer|Hex2Cer|GlcCer|TAG|TG|DAG|DG|MAG|MG|CE|ACar|CAR|FA)[ -]"
b <- rows[rows$file == "B", ]
lip <- b[grepl(lip_pat, trimws(b$original_name)), ]
message("Lipid rows in B: ", nrow(lip))
set.seed(20261007)
lp <- lip[sort(sample(nrow(lip), min(40, nrow(lip)))), ]
write.csv(data.frame(row_b = lp$row_id, name_b = lp$original_name, inchikey_b = lp$inchikey,
                     matches_in_A = "", relationship = "", evidence = "", stringsAsFactors = FALSE),
          file.path(out, "label_lipids.csv"), row.names = FALSE, na = "")
message("Lipid rows drawn: ", nrow(lp))
message("Sheets written to ", out)
