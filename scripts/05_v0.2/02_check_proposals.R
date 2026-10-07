# Runs propose_lipid_matches() through the full audit and compares its
# output with the proposals reviewed by hand in
# notes/v02_proposal_review.csv (made by scripts/05_v0.2/01_feasibility.R).
#
# Two things can differ:
#   - a reviewed proposal is missing: the package lost something the
#     feasibility logic found. This should not happen.
#   - a proposal is new: most likely from group-name rows, which get a
#     lipid profile since v0.2. New proposals have not been reviewed yet.
#
# Output: results/v02_check/
# Run from the project root: Rscript scripts/05_v0.2/02_check_proposals.R

pkgload::load_all(".", quiet = TRUE)

base <- "data-raw/MetLinkR/inst/extdata/Example_output"
files <- data.frame(
  path         = file.path(base, c("LEOCC_Metabolon_Annotations_metLinkR.csv",
                                   "Broad_2022Aug_annotations_metLinkR.csv")),
  name_col     = c("CHEMICAL_NAME", "Metabolite"),
  inchikey_col = c("INCHIKEY", "inChIKey"),
  label        = c("LEOCC_Metabolon", "Broad_2022Aug"),
  stringsAsFactors = FALSE
)

audit <- audit_metlinkr(files,
                        swisslipids  = "data-raw/swisslipids/lipids.tsv",
                        expected_md5 = "b6ddd5bc3003967dce95e905e5748953",
                        encoding     = "Latin-1")
m <- audit$lipid_matches

reviewed <- read.csv("notes/v02_proposal_review.csv", colClasses = "character")
key_now  <- paste(m$row_a, m$row_b)
key_rev  <- paste(reviewed$row_a, reviewed$row_b)

missing  <- reviewed[!key_rev %in% key_now, c("row_a", "row_b", "name_a", "name_b")]
new      <- m[!key_now %in% key_rev, c("row_a", "row_b", "name_a", "name_b",
                                       "lipid_species", "relationship", "rule_id",
                                       "shape", "already_grouped")]

out_dir <- "results/v02_check"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(m, file.path(out_dir, "lipid_matches.csv"), row.names = FALSE)
write.csv(new, file.path(out_dir, "new_proposals.csv"), row.names = FALSE)

message("Proposals from the package: ", nrow(m))
message("Reviewed proposals missing from the package output: ", nrow(missing))
if (nrow(missing) > 0) print(missing, row.names = FALSE)
message("Proposals not reviewed yet: ", nrow(new))
if (nrow(new) > 0) print(new[, c("name_a", "name_b", "relationship", "rule_id", "shape")],
                         row.names = FALSE)
message("\nLipid profile sources:")
print(table(audit$rows$lipid_source, useNA = "ifany"))
