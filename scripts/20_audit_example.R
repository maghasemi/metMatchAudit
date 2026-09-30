# Runs the v0.1 audit on the MetLinkR example output (LEOCC Metabolon and
# Broad 2022 files, the two with an InChIKey column).
#
# The manual review in notes/review.csv is not finished. To keep it
# blind, the groups listed there are held back:
#   results/audit_v0.1/          everything except the reviewed groups
#   results/audit_v0.1_sealed.rds  the full audit, for the evaluation
#                                  after the review is done. Do not open
#                                  it before then.
#
# Needs: scripts/00_get_data.sh and the SwissLipids file. Structure
# checks (R09b, C04) run when RDKit is available through reticulate.
# Run from the project root: Rscript scripts/20_audit_example.R

pkgload::load_all(".", quiet = TRUE)

base <- "data-raw/MetLinkR/inst/extdata/Example_output"
files <- data.frame(
  path         = file.path(base, c("LEOCC_Metabolon_Annotations_metLinkR.csv",
                                   "Broad_2022Aug_annotations_metLinkR.csv")),
  name_col     = c("CHEMICAL_NAME", "Metabolite"),
  inchikey_col = c("INCHIKEY", "inChIKey"),
  smiles_col   = c("SMILES", "smiles"),
  inchi_col    = c(NA, "inChI"),
  label        = c("LEOCC_Metabolon", "Broad_2022Aug"),
  stringsAsFactors = FALSE
)

audit <- audit_metlinkr(files,
                        swisslipids  = "data-raw/swisslipids/lipids.tsv",
                        expected_md5 = "b6ddd5bc3003967dce95e905e5748953",
                        encoding     = "Latin-1")

saveRDS(audit, "results/audit_v0.1_sealed.rds")

reviewed <- read.csv("notes/review.csv", colClasses = "character")$metlinkr_name
open <- audit
open$pairs  <- audit$pairs[!audit$pairs$group_id %in% reviewed, ]
open$groups <- audit$groups[!audit$groups$group_id %in% reviewed, ]
held_rows   <- audit$rows$row_id[audit$rows$group_id %in% reviewed]
open$rows   <- audit$rows[!audit$rows$row_id %in% held_rows, ]
open$row_flags <- audit$row_flags[!audit$row_flags$row_id %in% held_rows, ]

write_audit(open, "results/audit_v0.1")

message("Held back ", length(unique(reviewed)), " reviewed groups.")
print(open)
