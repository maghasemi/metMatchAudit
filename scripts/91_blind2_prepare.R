# Second blind test, step 2: CSV copies of the two annotation workbooks.
#
# Every column is read as text so nothing is converted on the way. The
# only change is one column name in B ("Pubchem CID" becomes PUBCHEM), so
# MetLinkR and this package read it without spaces. Plan:
# notes/decisions.md, 2026-10-07, second blind test.
#
# Input and output stay in data-raw/own/, which git ignores.
# Run from the project root: Rscript scripts/91_blind2_prepare.R

library(readxl)

dir <- "data-raw/own"

prepare <- function(src, out, expected_rows, rename = NULL) {
  d <- read_excel(file.path(dir, src), col_types = "text")
  if (nrow(d) != expected_rows) stop(src, " has ", nrow(d), " rows, expected ", expected_rows, ".")
  for (old in names(rename)) names(d)[names(d) == old] <- rename[[old]]
  write.csv(d, file.path(dir, out), row.names = FALSE, na = "")
  message(out, ": ", nrow(d), " rows, columns ", paste(names(d), collapse = ", "),
          ", md5 ", unname(tools::md5sum(file.path(dir, out))))
}

prepare("MyData.xlsx",  "A_metabolon.csv", 494)
prepare("2ndData.xlsx", "B_ding2021.csv",  1547, rename = c("Pubchem CID" = "PUBCHEM"))
