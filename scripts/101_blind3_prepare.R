# Third blind test, step 2: CSV copies of the two Lai et al. 2021
# workbooks (Supplementary Data 2 and 3). Plan: notes/decisions.md,
# 2026-10-07, third blind test.
#
# Every column is read as text. Only the compound name and the PubChem CID
# are kept, renamed BIOCHEMICAL and PUBCHEM as for the other files; MSI
# level is kept for reference. Nothing else changes.
#
# The workbooks are open access (CC BY 4.0) and live in data-raw/lai2021/
# with their CSV copies, all committed. The md5 of the private file A is
# printed so it can be checked against the second test.
#
# Run from the project root: Rscript scripts/101_blind3_prepare.R

library(readxl)

dir <- "data-raw/lai2021"
md5 <- function(p) unname(tools::md5sum(p))

prepare <- function(src, src_md5, name_col, out, expected_rows) {
  path <- file.path(dir, src)
  if (md5(path) != src_md5) stop(src, " has md5 ", md5(path), ", expected ", src_md5, ".")
  d <- read_excel(path, skip = 2, col_types = "text")
  d <- d[!is.na(d[[name_col]]), ]
  if (nrow(d) != expected_rows) stop(src, " has ", nrow(d), " rows, expected ", expected_rows, ".")
  out_d <- data.frame(BIOCHEMICAL = d[[name_col]], PUBCHEM = d[["PubChem_CID"]],
                      MSI_LEVEL = d[["MSI level"]], stringsAsFactors = FALSE)
  write.csv(out_d, file.path(dir, out), row.names = FALSE, na = "")
  message(out, ": ", nrow(out_d), " rows, ", sum(!is.na(out_d$PUBCHEM)), " with a CID, md5 ",
          md5(file.path(dir, out)))
}

prepare("41467_2021_26209_MOESM4_ESM.xlsx", "24b6e7dcf7e99befe9a49aa222dca921",
        "Compound name (feces)", "C1_lai_feces.csv", 533)
prepare("41467_2021_26209_MOESM5_ESM.xlsx", "22ca253e0e86ffdcb6fb6046f3794667",
        "Compound name (serum)", "C2_lai_serum.csv", 231)

a <- "data-raw/own/A_metabolon.csv"
message("A_metabolon.csv md5 ", md5(a), " (compare with the second test, scripts/91_blind2_prepare.R)")
