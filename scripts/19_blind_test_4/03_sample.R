# Fourth blind test, step 4: draws the groups and lipid rows to label and
# writes the two labeling sheets. Seed and rules are in notes/decisions.md
# (2026-11-08, fourth blind test, plan). Reads only MetLinkR's output and
# the names in it; the audit is not run here.
#
# The files are CC BY-NC-ND, so the sheets stay in data-raw/blind4/labels/,
# which git ignores. Their md5s go into notes/blind4_labels.md once the
# labels are accepted.
#
# Run from the project root: Rscript scripts/19_blind_test_4/03_sample.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.11.0")

mout <- "data-raw/blind4/metlinkr_output"
out  <- "data-raw/blind4/labels"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (file.exists(file.path(out, "label_pairs.csv"))) stop("Sheets exist already. Delete them on purpose to redraw.")

files <- data.frame(
  path         = file.path(mout, c("F1_metLinkR.csv", "F2_metLinkR.csv", "F3_metLinkR.csv")),
  name_col     = "BIOCHEMICAL",
  inchikey_col = NA,
  pubchem_col  = c("PUBCHEM", NA, NA),
  hmdb_col     = c("HMDB", NA, NA),
  label        = c("F1", "F2", "F3"),
  stringsAsFactors = FALSE
)
rows <- read_metlinkr(files)
print(table(rows$file))

# Groups with rows from at least two files, in a fixed order before sampling.
# A group is in the "same" stratum when all its names are one name after
# normalize_name(), and in the "different" stratum otherwise.
grp <- rows[!is.na(rows$group_id), ]
by_group <- split(grp, grp$group_id)
multi <- sort(names(by_group)[vapply(by_group, function(d) length(unique(d$file)) >= 2, logical(1))])
n_names <- vapply(by_group[multi], function(d) length(unique(normalize_name(d$original_name))), integer(1))
same <- multi[n_names == 1]
differ <- multi[n_names > 1]
combo <- vapply(by_group[multi], function(d) paste(sort(unique(d$file)), collapse = "+"), character(1))
message("MetLinkR groups with rows from at least two files: ", length(multi),
        " (same name: ", length(same), ", different names: ", length(differ), ")")
print(table(combo, stratum = ifelse(n_names == 1, "same", "different")))
if (length(multi) < 30) message("Fewer than 30 groups: the report must say the sample is too small for the thresholds.")

set.seed(20261108)
draw <- function(x, n) if (length(x) <= n) x else sample(x, n)
picked_same <- draw(same, 30)
picked_diff <- draw(differ, 30)
if (length(same) < 30) message("Same-name stratum has only ", length(same), " groups; all taken.")
if (length(differ) < 30) message("Different-name stratum has only ", length(differ), " groups; all taken.")
picked <- c(picked_same, picked_diff)
stratum <- setNames(rep(c("same", "different"), c(length(picked_same), length(picked_diff))), picked)

pairs <- make_pairs(rows[rows$group_id %in% picked, ])
r <- function(id, col) rows[[col]][match(id, rows$row_id)]
sheet <- data.frame(
  group_id = pairs$group_id,
  stratum  = unname(stratum[pairs$group_id]),
  row_a = pairs$row_a, name_a = r(pairs$row_a, "original_name"),
  pubchem_a = r(pairs$row_a, "pubchem"), hmdb_a = r(pairs$row_a, "hmdb"),
  row_b = pairs$row_b, name_b = r(pairs$row_b, "original_name"),
  pubchem_b = r(pairs$row_b, "pubchem"), hmdb_b = r(pairs$row_b, "hmdb"),
  cross_file = r(pairs$row_a, "file") != r(pairs$row_b, "file"),
  category = "", evidence = "",
  stringsAsFactors = FALSE
)
sheet <- sheet[order(sheet$stratum != "same", sheet$group_id, !sheet$cross_file), ]
write.csv(sheet, file.path(out, "label_pairs.csv"), row.names = FALSE, na = "")
message("Groups drawn: ", length(picked), " (same ", length(picked_same), ", different ",
        length(picked_diff), "); pairs to label: ", nrow(sheet),
        " (", sum(sheet$cross_file), " across files)")

# Lipid and carnitine rows of all three files, by the pattern of the third test.
lip_pat <- "^((Lyso)?(PC|PE|PS|PI|PG|PA)|SM|Cer|CerP|GlcCer|GalCer|LacCer|TG|DG|MG|CE)\\(|carnitine$"
cand <- rows[grepl(lip_pat, trimws(rows$original_name), ignore.case = TRUE, perl = TRUE), ]
cand <- cand[order(cand$row_id), ]
lip <- cand[sort(draw(seq_len(nrow(cand)), 40)), ]
message("Lipid and carnitine rows matching the pattern: ", nrow(cand), "; taken: ", nrow(lip))
print(table(lip$file))
write.csv(data.frame(row = lip$row_id, name = lip$original_name, pubchem = lip$pubchem,
                     hmdb = lip$hmdb, matches_in_other_files = "", relationship = "",
                     evidence = "", stringsAsFactors = FALSE),
          file.path(out, "label_lipids.csv"), row.names = FALSE, na = "")
message("Sheets written to ", out)
