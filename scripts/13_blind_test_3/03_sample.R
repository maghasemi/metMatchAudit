# Third blind test, step 4: draws the groups and lipid rows to label and
# writes the two labeling sheets. Seed and rules are in notes/decisions.md
# (2026-10-07, third blind test). Reads only MetLinkR's output; the audit
# is not run here.
#
# The sheets hold names from the private file A, so they stay in
# data-raw/own/blind3/. Their md5 goes into notes/blind3_labels.md once
# they are filled in.
#
# Run from the project root: Rscript scripts/13_blind_test_3/03_sample.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.7.0")

mout <- "data-raw/own/metlinkr_output_blind3"
out  <- "data-raw/own/blind3"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (file.exists(file.path(out, "label_pairs.csv"))) stop("Sheets exist already. Delete them on purpose to redraw.")

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
rows <- read_metlinkr(files)
print(table(rows$file))

# Groups with rows from at least two files, in a fixed order before sampling.
g <- split(rows$file, rows$group_id)
multi <- sort(names(g)[vapply(g, function(f) length(unique(f)) >= 2, logical(1))])
combo <- vapply(g[multi], function(f) paste(sort(unique(f)), collapse = "+"), character(1))
message("MetLinkR groups with rows from at least two files: ", length(multi))
print(table(combo))
if (length(multi) < 30) message("Fewer than 30 groups: the report must say the sample is too small for the thresholds.")

set.seed(20261008)
picked <- if (length(multi) <= 60) multi else sample(multi, 60)

pairs <- make_pairs(rows[rows$group_id %in% picked, ])
r <- function(id, col) rows[[col]][match(id, rows$row_id)]
sheet <- data.frame(
  group_id = pairs$group_id,
  row_a = pairs$row_a, name_a = r(pairs$row_a, "original_name"),
  pubchem_a = r(pairs$row_a, "pubchem"), hmdb_a = r(pairs$row_a, "hmdb"),
  row_b = pairs$row_b, name_b = r(pairs$row_b, "original_name"),
  pubchem_b = r(pairs$row_b, "pubchem"), hmdb_b = r(pairs$row_b, "hmdb"),
  cross_file = r(pairs$row_a, "file") != r(pairs$row_b, "file"),
  category = "", evidence = "",
  stringsAsFactors = FALSE
)
sheet <- sheet[order(sheet$group_id, !sheet$cross_file), ]
write.csv(sheet, file.path(out, "label_pairs.csv"), row.names = FALSE, na = "")
message("Groups drawn: ", length(picked), "; pairs to label: ", nrow(sheet),
        " (", sum(sheet$cross_file), " across files)")

# Lipid and carnitine rows of C1 and C2, by the pattern fixed in the plan.
lip_pat <- "^((Lyso)?(PC|PE|PS|PI|PG|PA)|SM|Cer|CerP|GlcCer|GalCer|LacCer|TG|DG|MG|CE)\\(|carnitine$"
cc <- rows[rows$file %in% c("C1", "C2"), ]
lip <- cc[grepl(lip_pat, trimws(cc$original_name), ignore.case = TRUE, perl = TRUE), ]
message("Lipid and carnitine rows in C1 and C2: ", nrow(lip), " (plan: 40, all taken)")
write.csv(data.frame(row = lip$row_id, name = lip$original_name, pubchem = lip$pubchem,
                     matches_in_other_files = "", relationship = "", evidence = "",
                     stringsAsFactors = FALSE),
          file.path(out, "label_lipids.csv"), row.names = FALSE, na = "")
message("Sheets written to ", out)
