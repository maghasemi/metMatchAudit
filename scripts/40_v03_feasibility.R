# Feasibility test for v0.3: can PubChem turn a SMILES into InChI and
# InChIKey reliably enough to use for stereo comparison?
#
# The reference is scripts/40_rdkit_reference.csv: InChIKeys computed
# from the same SMILES with RDKit (scripts/40_rdkit_reference.py). It
# holds five named rows from the manual review, four of them structures
# that PubChem did not have a record for, and a random sample of 60 rows.
#
# Criteria are in notes/decisions.md (2026-10-01, v0.3 feasibility),
# written before this script was first run.
#
# Also counts how many rows have a structure (SMILES or InChI) and how
# many R11 pairs have one on both sides.
#
# Needs: install.packages(c("httr2", "jsonlite")), internet access.
# Output: results/v03_feasibility/
# Run from the project root: Rscript scripts/40_v03_feasibility.R

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
})

out_dir <- "results/v03_feasibility"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache_path <- file.path(out_dir, "pubchem_cache.rds")


# ---- PubChem -----------------------------------------------------------------
# One POST request per SMILES, so characters such as "/" and "#" in the
# SMILES do not need escaping. Answers are cached by SMILES, so a rerun
# does not ask PubChem again.

pug_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/InChI,InChIKey/JSON"

query_pubchem <- function(smiles) {
  started <- Sys.time()
  resp <- tryCatch(
    request(pug_url) |>
      req_body_form(smiles = smiles) |>
      req_user_agent("metMatchAudit feasibility test") |>
      req_timeout(30) |>
      req_error(is_error = function(r) FALSE) |>
      req_perform(),
    error = function(e) e
  )
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  empty <- list(status = NA_integer_, cid = NA_character_, inchi = NA_character_,
                inchikey = NA_character_, message = NA_character_, seconds = seconds)
  if (inherits(resp, "error")) {
    empty$message <- conditionMessage(resp)
    return(empty)
  }
  empty$status <- resp_status(resp)
  body <- tryCatch(fromJSON(resp_body_string(resp)), error = function(e) NULL)
  props <- body$PropertyTable$Properties
  if (is.null(props) || nrow(props) == 0) {
    empty$message <- if (!is.null(body$Fault$Message)) body$Fault$Message else "no properties returned"
    return(empty)
  }
  empty$cid      <- if ("CID" %in% names(props)) as.character(props$CID[1]) else NA_character_
  empty$inchi    <- if ("InChI" %in% names(props)) props$InChI[1] else NA_character_
  empty$inchikey <- if ("InChIKey" %in% names(props)) props$InChIKey[1] else NA_character_
  empty$n_records <- nrow(props)
  empty
}

cache <- if (file.exists(cache_path)) readRDS(cache_path) else list()

ref <- read.csv("scripts/40_rdkit_reference.csv", colClasses = "character")
# Drop the CXSMILES extension, as the RDKit script does.
ref$smiles_clean <- sub(" .*$", "", ref$smiles)

for (s in unique(ref$smiles_clean)) {
  if (!is.null(cache[[s]])) next
  cache[[s]] <- query_pubchem(s)
  saveRDS(cache, cache_path)
  Sys.sleep(0.3)  # PubChem asks for no more than 5 requests per second
}

get_field <- function(field) {
  vapply(ref$smiles_clean, function(s) {
    v <- cache[[s]][[field]]
    if (is.null(v)) NA_character_ else as.character(v)
  }, character(1), USE.NAMES = FALSE)
}
ref$pubchem_status  <- get_field("status")
ref$pubchem_cid     <- get_field("cid")
ref$pubchem_key     <- get_field("inchikey")
ref$pubchem_inchi   <- get_field("inchi")
ref$pubchem_message <- get_field("message")
ref$seconds         <- as.numeric(get_field("seconds"))

ref$match_25   <- !is.na(ref$pubchem_key) & substr(ref$pubchem_key, 1, 25) == substr(ref$rdkit_key, 1, 25)
ref$match_full <- !is.na(ref$pubchem_key) & ref$pubchem_key == ref$rdkit_key

write.csv(ref, file.path(out_dir, "pubchem_vs_rdkit.csv"), row.names = FALSE)


# ---- Coverage -------------------------------------------------------------------

base <- "data-raw/MetLinkR/inst/extdata/Example_output"
leocc <- read.csv(file.path(base, "LEOCC_Metabolon_Annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))
broad <- read.csv(file.path(base, "Broad_2022Aug_annotations_metLinkR.csv"),
                  colClasses = "character", fileEncoding = "latin1", na.strings = c("", "NA"))
structure_rows <- rbind(
  data.frame(row_id = paste0("LEOCC_Metabolon:", seq_len(nrow(leocc))),
             has_smiles = !is.na(leocc$SMILES), has_inchi = FALSE),
  data.frame(row_id = paste0("Broad_2022Aug:", seq_len(nrow(broad))),
             has_smiles = !is.na(broad$smiles), has_inchi = !is.na(broad$inChI))
)
structure_rows$has_structure <- structure_rows$has_smiles | structure_rows$has_inchi

audit <- readRDS("results/audit_v0.1_sealed.rds")
r11 <- audit$pairs[audit$pairs$rule_id == "R11", ]
has <- setNames(structure_rows$has_structure, structure_rows$row_id)
r11$both_have_structure <- has[r11$row_a] & has[r11$row_b]
write.csv(r11, file.path(out_dir, "r11_pairs.csv"), row.names = FALSE)


# ---- Summary -----------------------------------------------------------------------

named <- ref[ref$set == "named", ]
message("Named rows (criterion: 5 of 5 agree on the first 25 characters):")
print(named[, c("name", "rdkit_key", "pubchem_key", "pubchem_cid", "match_25", "pubchem_message")],
      row.names = FALSE)

message("\nSample rows agreeing on the first 25 characters: ",
        sum(ref$match_25[ref$set == "sample"]), " of ", sum(ref$set == "sample"))
message("Rows with no answer from PubChem: ", sum(is.na(ref$pubchem_key)))
bad <- ref[!ref$match_25, c("set", "name", "rdkit_key", "pubchem_key", "pubchem_message")]
if (nrow(bad) > 0) {
  message("Rows that do not agree:")
  print(bad, row.names = FALSE)
}

message("\nSeconds per request: median ", round(median(ref$seconds, na.rm = TRUE), 2),
        ", max ", round(max(ref$seconds, na.rm = TRUE), 2))
message("\nRows with a structure: LEOCC ",
        sum(structure_rows$has_structure[grepl("^LEOCC", structure_rows$row_id)]), " of ", nrow(leocc),
        ", Broad ", sum(structure_rows$has_structure[grepl("^Broad", structure_rows$row_id)]), " of ", nrow(broad))
message("R11 pairs with a structure on both sides: ", sum(r11$both_have_structure), " of ", nrow(r11))
message("\nFiles written to ", out_dir)
