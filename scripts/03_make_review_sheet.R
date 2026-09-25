# Creates notes/review.csv, the sheet you fill in by hand during review.
#
# Which groups go in:
#   every "different_connectivity" group
#   every same-connectivity group ("same_connectivity_different_block2"
#   or "protonation_only") if there are 40 or fewer; otherwise a random
#   sample of 40 with a fixed seed, so the same sample comes out each run
#
# The script refuses to overwrite an existing review.csv, because that
# file holds your manual work.
#
# Fill in two columns per line:
#   category      one of A to F, as defined in notes/decisions.md
#   evidence_note what the PubChem records showed, in a few words
# Run from the project root: Rscript scripts/03_make_review_sheet.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

out_path <- "notes/review.csv"

if (file.exists(out_path)) {
  stop(out_path, " already exists. Rename or delete it by hand if you ",
       "really want a new one.")
}

groups  <- read_csv("results/inchikey_groups.csv",
                    col_types = cols(.default = col_character()),
                    show_col_types = FALSE)
flagged <- read_csv("results/flagged_rows.csv",
                    col_types = cols(.default = col_character()),
                    show_col_types = FALSE)

diff_conn <- filter(groups, pattern == "different_connectivity")
same_conn <- filter(groups, pattern %in% c("same_connectivity_different_block2",
                                           "protonation_only"))

sample_size <- 40
if (nrow(same_conn) > sample_size) {
  set.seed(20260923)
  same_conn <- slice_sample(same_conn, n = sample_size)
}

selected <- bind_rows(diff_conn, same_conn)

# One line per group. The InChIKey search link opens the PubChem
# record list for that key.
review <- flagged |>
  semi_join(selected, by = "metlinkr_name") |>
  group_by(metlinkr_name, pattern) |>
  summarise(row_ids        = paste(row_id, collapse = " | "),
            original_names = paste(original_name, collapse = " | "),
            inchikeys      = paste(unique(inchikey), collapse = " | "),
            pubchem_links  = paste0("https://pubchem.ncbi.nlm.nih.gov/#query=",
                                    unique(inchikey), collapse = " | "),
            .groups        = "drop") |>
  mutate(category = "", evidence_note = "") |>
  arrange(pattern, metlinkr_name)

write_csv(review, out_path, na = "")
message("Wrote ", out_path, " with ", nrow(review), " groups to review")
