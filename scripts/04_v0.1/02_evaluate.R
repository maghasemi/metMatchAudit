# Compares the v0.1 audit with the manual review, group by group.
#
# Run this only after notes/review.csv is complete and committed. It is
# the first script allowed to open results/audit_v0.1_sealed.rds.
#
# The review gives one category per group, for the most severe problem
# (E > C > D > B > A). The audit gives one relationship per group, for
# its most severe pair. Each category maps to the relationship the tool
# should report:
#
#   A exact   B stereo_variant   C distinct   D broader_narrower
#   F unresolved
#   E has no expected relationship in v0.1: source errors need the
#     checks planned for v0.3. For E groups we only ask whether the tool
#     called them exact, which would endorse a merge the review rejected.
#
# Every group gets one outcome:
#   agree       the tool reports the expected relationship
#   safe_miss   the tool says unresolved where the review has an answer
#   wrong       the tool claims a different definite relationship
#   E_not_exact / E_exact   for E groups, see above
#
# Output: results/evaluation_v0.1/
# Run from the project root: Rscript scripts/04_v0.1/02_evaluate.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

audit  <- readRDS("results/audit_v0.1_sealed.rds")
review <- read_csv("notes/review.csv", col_types = cols(.default = "c"))

stopifnot(nrow(review) == 53, all(review$category %in% LETTERS[1:6]))

expected_map <- c(A = "exact", B = "stereo_variant", C = "distinct",
                  D = "broader_narrower", E = NA, F = "unresolved")

tool <- audit$groups |>
  select(group_id, tool_relationship = worst_relationship, tool_rule = worst_rule,
         n_rows, n_pairs)

eval <- review |>
  select(metlinkr_name, category, evidence_note) |>
  left_join(tool, by = c("metlinkr_name" = "group_id")) |>
  mutate(
    expected = unname(expected_map[category]),
    outcome = case_when(
      is.na(tool_relationship)                          ~ "not_in_audit",
      category == "E" & tool_relationship == "exact"    ~ "E_exact",
      category == "E"                                   ~ "E_not_exact",
      tool_relationship == expected                     ~ "agree",
      tool_relationship == "unresolved"                 ~ "safe_miss",
      TRUE                                              ~ "wrong"
    )
  )

# Every reviewed group should be in the audit, since both come from the
# same two files. Stop if not, rather than report a partial comparison.
if (any(eval$outcome == "not_in_audit")) {
  print(filter(eval, outcome == "not_in_audit") |> select(metlinkr_name))
  stop("Some reviewed groups are missing from the audit.")
}

by_category <- eval |> count(category, outcome, name = "groups")
overall     <- eval |> count(outcome, name = "groups")
by_rule     <- eval |> count(category, tool_rule, name = "groups")

disagreements <- eval |>
  filter(outcome != "agree") |>
  arrange(outcome, category, metlinkr_name) |>
  select(metlinkr_name, category, expected, tool_relationship, tool_rule,
         outcome, evidence_note)

out_dir <- "results/evaluation_v0.1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(eval,          file.path(out_dir, "per_group.csv"))
write_csv(by_category,   file.path(out_dir, "by_category.csv"))
write_csv(overall,       file.path(out_dir, "overall.csv"))
write_csv(by_rule,       file.path(out_dir, "by_rule.csv"))
write_csv(disagreements, file.path(out_dir, "disagreements.csv"))

message("Overall:")
print(overall)
message("\nBy review category:")
print(by_category, n = Inf)
message("\nRule behind the tool's answer, per review category:")
print(by_rule, n = Inf)
message("\nFiles written to ", out_dir)
