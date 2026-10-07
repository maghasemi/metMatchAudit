# Reads the SwissLipids download once and saves the columns we need
# as a smaller R file, so later scripts do not re-read 690 MB of text.
#
# Before running: download lipids.tsv from
# https://www.swisslipids.org/#/downloads and put it in
# data-raw/swisslipids/. The md5 check below stops the script if the
# file differs from the one this project was built on. If you download
# a newer release on purpose, log it in notes/decisions.md and update
# expected_md5.
#
# Output:
#   data-raw/swisslipids/lipids_slim.rds
#   results/lipid_feasibility/swisslipids_levels.csv
# Run from the project root: Rscript scripts/01_setup/02_swisslipids_slim.R

suppressPackageStartupMessages(library(data.table))

path         <- "data-raw/swisslipids/lipids.tsv"
out_path     <- "data-raw/swisslipids/lipids_slim.rds"
out_dir      <- "results/lipid_feasibility"
expected_md5 <- "b6ddd5bc3003967dce95e905e5748953"

if (!file.exists(path)) stop("Cannot find ", path, ". See the note at the top of this script.")

md5 <- unname(tools::md5sum(path))
if (md5 != expected_md5) {
  stop("md5 of ", path, " is ", md5, ", expected ", expected_md5, ".")
}

# Check the header by position. Column names carry spaces, asterisks
# and in two cases a trailing space, so we compare trimmed names and
# then select by position.
header <- strsplit(readLines(path, n = 1), "\t", fixed = TRUE)[[1]]
wanted <- c("1"  = "Lipid ID",
            "2"  = "Level",
            "3"  = "Name",
            "4"  = "Abbreviation*",
            "6"  = "Lipid class*",
            "7"  = "Parent",
            "11" = "InChI key (pH7.3)",
            "12" = "Formula (pH7.3)",
            "25" = "CHEBI",
            "26" = "LIPID MAPS",
            "27" = "HMDB")
pos <- as.integer(names(wanted))
found <- trimws(header[pos])
if (!identical(found, unname(wanted))) {
  print(data.frame(position = pos, expected = unname(wanted), found = found))
  stop("SwissLipids header does not match what this script expects.")
}

# quote = "" because names and SMILES can contain quote characters that
# are not meant as CSV quoting.
sl <- fread(path, sep = "\t", quote = "", select = pos,
            colClasses = "character", na.strings = "",
            encoding = "UTF-8", showProgress = FALSE)
setnames(sl, c("slm_id", "level", "name", "abbrev", "lipid_class", "parent",
               "inchikey_raw", "formula", "chebi", "lipidmaps", "hmdb"))

# Keys come as "InChIKey=XXXX-YYYY-Z". Strip the prefix and keep only
# values with the standard format.
sl[, inchikey := sub("^InChIKey=", "", inchikey_raw)]
sl[!grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", inchikey), inchikey := NA_character_]
sl[, inchikey_raw := NULL]
sl[, key25 := substr(inchikey, 1, 25)]
sl[, key14 := substr(inchikey, 1, 14)]

# Some rows might list more than one parent. The feasibility script
# handles that, but we want to know how common it is.
multi_parent <- sl[grepl("[|,]", parent), .N]

levels_table <- sl[, .(n = .N, with_inchikey = sum(!is.na(inchikey)),
                       with_abbrev = sum(!is.na(abbrev))), by = level][order(-n)]

message("Rows read: ", nrow(sl))
message("Rows with a valid InChIKey: ", sl[!is.na(inchikey), .N])
message("Rows with more than one parent listed: ", multi_parent)
print(levels_table)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(levels_table, file.path(out_dir, "swisslipids_levels.csv"))
saveRDS(sl, out_path)
message("Wrote ", out_path)
