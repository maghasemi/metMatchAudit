# Builds one row-level table from two MetLinkR example output files.
#
# We only use the two files that carry an InChIKey column:
#   LEOCC_Metabolon_Annotations_metLinkR.csv  (Metabolon, column INCHIKEY)
#   Broad_2022Aug_annotations_metLinkR.csv    (Broad, column inChIKey)
# Each file is the original input returned by MetLinkR with one extra
# column, "Standardized name", holding the name MetLinkR assigned.
#
# Output: results/rows.csv, one line per input row.
# Run from the project root: Rscript scripts/01_build_rows.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

base <- "data-raw/MetLinkR/inst/extdata/Example_output"

if (!dir.exists(base)) {
  stop("Cannot find ", base, ". Run scripts/00_get_data.sh first.")
}

# Read everything as text. Numeric parsing would mangle IDs such as
# PubChem CIDs with leading zeros or values like "6950182;99288".
read_chr <- function(path) {
  read_csv(path,
           col_types = cols(.default = col_character()),
           locale = locale(encoding = "latin1"),
           na = c("", "NA"),
           show_col_types = FALSE)
}

leocc <- read_chr(file.path(base, "LEOCC_Metabolon_Annotations_metLinkR.csv"))
broad <- read_chr(file.path(base, "Broad_2022Aug_annotations_metLinkR.csv"))

# These counts match the files at commit a6aa145. If they differ, the
# download is wrong or the files changed upstream. Stop and check.
stopifnot(nrow(leocc) == 1548, nrow(broad) == 511)

rows <- bind_rows(
  leocc |> transmute(file          = "LEOCC_Metabolon",
                     row_in_file   = row_number(),
                     original_name = CHEMICAL_NAME,
                     inchikey      = INCHIKEY,
                     metlinkr_name = `Standardized name`),
  broad |> transmute(file          = "Broad_2022Aug",
                     row_in_file   = row_number(),
                     original_name = Metabolite,
                     inchikey      = inChIKey,
                     metlinkr_name = `Standardized name`)
) |>
  mutate(row_id   = paste(file, row_in_file, sep = ":"),
         inchikey = str_trim(inchikey),
         # A standard InChIKey is 14 letters, a hyphen, 10 letters,
         # a hyphen, and 1 letter. Anything else gets flagged here and
         # left as it is. We never edit source values.
         ik_valid = !is.na(inchikey) &
                    str_detect(inchikey, "^[A-Z]{14}-[A-Z]{10}-[A-Z]$")) |>
  relocate(row_id)

summary_table <- rows |>
  group_by(file) |>
  summarise(n_rows      = n(),
            has_ik      = sum(!is.na(inchikey)),
            invalid_ik  = sum(!is.na(inchikey) & !ik_valid),
            no_metlinkr = sum(is.na(metlinkr_name)),
            .groups     = "drop")

print(summary_table)

write_csv(rows, "results/rows.csv")
write_csv(summary_table, "results/rows_summary.csv")
message("Wrote results/rows.csv (", nrow(rows), " rows)")
