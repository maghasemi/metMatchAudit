# Feasibility test for v0.3 with RDKit through reticulate.
#
# 1. Bridge: RDKit called from R must give exactly the reference keys
#    (scripts/40_rdkit_reference.csv).
# 2. Stereo comparison: compare_stereo() (scripts/41_inchi_stereo_functions.R)
#    is run on every R11 pair of the audit. Where the group was reviewed,
#    the result is set against the manual label.
# 3. Row consistency (the planned check C04): for every row with a SMILES
#    or InChI, does the key computed from it match the key stored in the
#    same row?
#
# Criteria are in notes/decisions.md (2026-10-01, v0.3 RDKit
# feasibility), written before this script was first run.
#
# Needs: reticulate (>= 1.41), internet on the first run to install RDKit.
# Output: results/v03_rdkit/
# Run from the project root: Rscript scripts/41_v03_rdkit_feasibility.R

library(reticulate)
py_require("rdkit")
Chem  <- import("rdkit.Chem")
Inchi <- import("rdkit.Chem.inchi")
RDLogger <- import("rdkit.RDLogger")
RDLogger$DisableLog("rdApp.*")
source("scripts/41_inchi_stereo_functions.R")

out_dir <- "results/v03_rdkit"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# InChI from a SMILES. The CXSMILES extension after a space is dropped.
smiles_to_inchi <- function(smiles) {
  vapply(smiles, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    mol <- Chem$MolFromSmiles(sub(" .*$", "", s))
    if (is.null(mol)) return(NA_character_)
    out <- Inchi$MolToInchi(mol)
    if (is.null(out) || !nzchar(out)) NA_character_ else out
  }, character(1), USE.NAMES = FALSE)
}
inchi_to_key <- function(inchi) {
  vapply(inchi, function(x) {
    if (is.na(x)) return(NA_character_)
    k <- Inchi$InchiToInchiKey(x)
    if (is.null(k)) NA_character_ else k
  }, character(1), USE.NAMES = FALSE)
}


# ---- 1. Bridge ----------------------------------------------------------------

ref <- read.csv("scripts/40_rdkit_reference.csv", colClasses = "character")
ref$r_key <- inchi_to_key(smiles_to_inchi(ref$smiles))
bridge_ok <- sum(ref$r_key == ref$rdkit_key, na.rm = TRUE)


# ---- Structures for every row -------------------------------------------------
# Broad has its own InChI column, which is used as given. LEOCC has only
# SMILES, turned into InChI here. The structure is used only if its key
# agrees with the row's stored key on the first 25 characters; otherwise
# the row is a C04 case and its structure is not trusted.

base <- "data-raw/MetLinkR/inst/extdata/Example_output"
leocc <- read.csv(file.path(base, "LEOCC_Metabolon_Annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))
broad <- read.csv(file.path(base, "Broad_2022Aug_annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))

started <- Sys.time()
rows <- rbind(
  data.frame(row_id = paste0("LEOCC_Metabolon:", seq_len(nrow(leocc))),
             name = leocc$CHEMICAL_NAME, stored_key = leocc$INCHIKEY,
             smiles = leocc$SMILES, inchi = smiles_to_inchi(leocc$SMILES),
             inchi_source = ifelse(is.na(leocc$SMILES), NA, "rdkit_from_smiles")),
  data.frame(row_id = paste0("Broad_2022Aug:", seq_len(nrow(broad))),
             name = broad$Metabolite, stored_key = broad$inChIKey,
             smiles = broad$smiles, inchi = broad$inChI,
             inchi_source = ifelse(is.na(broad$inChI), NA, "file"))
)
rows$structure_key <- inchi_to_key(rows$inchi)
seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

# Also the key from the Broad SMILES, so the Broad SMILES can be checked
# against its own row too (xylose and quinic acid were found that way).
is_broad <- grepl("^Broad", rows$row_id)
rows$smiles_key <- NA_character_
rows$smiles_key[is_broad] <- inchi_to_key(smiles_to_inchi(rows$smiles[is_broad]))
rows$smiles_key[!is_broad] <- rows$structure_key[!is_broad]

k25 <- function(x) substr(x, 1, 25)
has_stored <- !is.na(rows$stored_key)
rows$c04_structure <- has_stored & !is.na(rows$structure_key) & k25(rows$structure_key) != k25(rows$stored_key)
rows$c04_smiles    <- has_stored & !is.na(rows$smiles_key) & k25(rows$smiles_key) != k25(rows$stored_key)
rows$usable_inchi  <- ifelse(!is.na(rows$structure_key) & !rows$c04_structure, rows$inchi, NA_character_)
write.csv(rows, file.path(out_dir, "row_structures.csv"), row.names = FALSE)


# ---- 2. Stereo comparison on R11 pairs ------------------------------------------

audit <- readRDS("results/audit_v0.1_sealed.rds")
r11 <- audit$pairs[audit$pairs$rule_id == "R11", c("group_id", "row_a", "row_b")]
inchi_of <- setNames(rows$usable_inchi, rows$row_id)
name_of  <- setNames(rows$name, rows$row_id)
r11$name_a <- name_of[r11$row_a]
r11$name_b <- name_of[r11$row_b]
r11$stereo <- mapply(function(a, b) {
  if (is.na(a) || is.na(b)) return("no_usable_structure")
  compare_stereo(a, b)
}, inchi_of[r11$row_a], inchi_of[r11$row_b])

review <- read.csv("notes/review.csv", colClasses = "character")
r11$review_category <- review$category[match(r11$group_id, review$metlinkr_name)]
write.csv(r11, file.path(out_dir, "r11_stereo.csv"), row.names = FALSE)


# ---- Summary -----------------------------------------------------------------------

message("1. Bridge: ", bridge_ok, " of ", nrow(ref), " reference keys reproduced exactly")
message("\n2. Stereo comparison on R11 pairs:")
print(table(stereo = r11$stereo, review = ifelse(is.na(r11$review_category), "not reviewed", r11$review_category)))
message("\nR11 pairs, one line each:")
print(r11[, c("review_category", "stereo", "name_a", "name_b")], row.names = FALSE)
rows$c04_block1 <- rows$c04_structure & substr(rows$structure_key, 1, 14) != substr(rows$stored_key, 1, 14)
message("\n3. Rows whose structure key differs from the stored key (first 25 characters): ",
        sum(rows$c04_structure), ", of which ", sum(rows$c04_block1), " differ in connectivity")
print(rows[rows$c04_structure, c("row_id", "name", "stored_key", "structure_key", "c04_block1")],
      row.names = FALSE)
# Broad SMILES often mark stereo as relative (CXSMILES "&1"). That part is
# dropped before RDKit reads them, so their stereo block cannot be
# compared with the stored key. Only the connectivity block is.
smiles_b1 <- is_broad & has_stored & !is.na(rows$smiles_key) &
  substr(rows$smiles_key, 1, 14) != substr(rows$stored_key, 1, 14)
message("Broad rows whose SMILES has a different connectivity block than the stored key: ",
        sum(smiles_b1))
print(rows[smiles_b1, c("row_id", "name", "stored_key", "smiles_key")], row.names = FALSE)
message("\nStructures computed in ", round(seconds, 1), " seconds")
message("Files written to ", out_dir)
