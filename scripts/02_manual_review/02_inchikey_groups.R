# Checks whether rows that MetLinkR gave the same standardized name
# also share the same structure, using the InChIKeys already in the files.
#
# An InChIKey has three parts, for example BPYKTIZUTYGOLE-KDUUSRDASA-N:
#   block1 (14 letters): the connectivity skeleton
#   block2 (10 letters): stereo, isotopes, and version flags
#   last letter:         protonation state
#
# The pattern column below only describes how the InChIKeys differ.
# It is not a relationship label yet. A block2 difference often means
# stereochemistry, but we confirm that during manual review, not here.
#
# Output:
#   results/inchikey_groups.csv   one line per MetLinkR name with 2+ rows
#   results/flagged_rows.csv      the rows behind every non-identical group
# Run from the project root: Rscript scripts/02_manual_review/02_inchikey_groups.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

rows <- read_csv("results/rows.csv",
                 col_types = cols(.default = col_character()),
                 show_col_types = FALSE)

ik <- rows |>
  filter(!is.na(metlinkr_name), ik_valid == "TRUE") |>
  mutate(block1 = str_sub(inchikey, 1, 14),
         block2 = str_sub(inchikey, 16, 25),
         prot   = str_sub(inchikey, 27, 27))

groups <- ik |>
  group_by(metlinkr_name) |>
  summarise(n_rows   = n(),
            n_files  = n_distinct(file),
            n_ik     = n_distinct(inchikey),
            n_block1 = n_distinct(block1),
            n_b1b2   = n_distinct(paste(block1, block2)),
            files    = paste(sort(unique(file)), collapse = ","),
            .groups  = "drop") |>
  filter(n_rows > 1) |>
  mutate(pattern = case_when(
    n_block1 > 1 ~ "different_connectivity",
    n_b1b2   > 1 ~ "same_connectivity_different_block2",
    n_ik     > 1 ~ "protonation_only",
    TRUE         ~ "identical_inchikey"
  )) |>
  arrange(pattern, metlinkr_name)

pattern_counts <- count(groups, pattern)
print(pattern_counts)

flagged <- ik |>
  semi_join(filter(groups, pattern != "identical_inchikey"),
            by = "metlinkr_name") |>
  left_join(select(groups, metlinkr_name, pattern), by = "metlinkr_name") |>
  arrange(pattern, metlinkr_name, file)

write_csv(groups, "results/inchikey_groups.csv")
write_csv(pattern_counts, "results/pattern_counts.csv")
write_csv(flagged, "results/flagged_rows.csv")
message("Groups with 2+ rows: ", nrow(groups))
message("Flagged rows written: ", nrow(flagged))
