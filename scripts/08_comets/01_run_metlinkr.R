# Runs MetLinkR on the five COMETS files that have expert mappings
# (ManualMappings.csv, MetLinkR SI), with MetLinkR's default settings.
# The output is the "MetLinkR alone" side of the COMETS external test of
# metMatchAudit, and the input that metMatchAudit audits.
#
# This project is kept apart from metMatchAudit so that MetLinkR, RaMP
# and Java stay out of the package's renv.lock.
#
# Run from the root of this project: Rscript run_metlinkr_comets.R
# Output: comets_run/metLinkR_output/ and comets_run/run_log.txt, then
# copied to metMatchAudit/data-raw/comets/metlinkr_output/.

suppressPackageStartupMessages({
  library(metLinkR)
  library(RaMP)
})

src     <- path.expand("~/Projects/metMatchAudit/data-raw/comets")
dest    <- file.path(src, "metlinkr_output")
run_dir <- "comets_run"

files <- c("inputs_1539353822.csv", "inputs_1539354673.csv", "inputs_1553799316.csv",
           "inputs_1555609921.csv", "inputs_1571898296.csv")
missing <- files[!file.exists(file.path(src, files))]
if (length(missing) > 0) stop("Cannot find in ", src, ": ", paste(missing, collapse = ", "))

dir.create(run_dir, showWarnings = FALSE)
file.copy(file.path(src, files), run_dir, overwrite = TRUE)

# Columns follow the headers of each file. Only two files carry an
# identifier: inputs_1539354673 has HMDB IDs, inputs_1555609921 has
# PubChem CIDs. All five have metabolite_name.
manifest <- data.frame(
  FileNames       = files,
  ShortFileName   = sub("\\.csv$", "", files),
  HMDB            = c(NA, "HMDB_ID", NA, NA, NA),
  Metabolite_Name = "metabolite_name",
  PubChem_CID     = c(NA, NA, NA, "PUBCHEM", NA),
  KEGG            = NA,
  LIPIDMAPS       = NA,
  chebi           = NA
)
write.csv(manifest, file.path(run_dir, "manifest.csv"), row.names = FALSE, na = "NA")

# MetLinkR reads the manifest and writes metLinkR_output/ relative to the
# working directory.
old <- setwd(run_dir)
on.exit(setwd(old), add = TRUE)

started <- Sys.time()
result <- harmonizeInputSheets(inputcsv = "manifest.csv", n_cores = 1)
minutes <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)

desc <- utils::packageDescription("metLinkR")
log <- c(
  paste("Run date:", format(started, "%Y-%m-%d %H:%M %Z")),
  paste("Run time (minutes):", minutes),
  paste("metLinkR version:", desc$Version),
  paste("metLinkR GitHub commit:", if (is.null(desc$RemoteSha)) "unknown" else desc$RemoteSha),
  paste("RaMP package version:", as.character(utils::packageVersion("RaMP"))),
  "RaMP database versions on this machine:",
  utils::capture.output(RaMP::listAvailableRaMPDbVersions()),
  "Settings: MetLinkR defaults, n_cores = 1",
  "",
  utils::capture.output(sessionInfo())
)
writeLines(log, "run_log.txt")

setwd(old)
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
file.copy(list.files(file.path(run_dir, "metLinkR_output"), full.names = TRUE), dest, overwrite = TRUE)
file.copy(file.path(run_dir, c("manifest.csv", "run_log.txt")), dest, overwrite = TRUE)
message("MetLinkR finished in ", minutes, " minutes. Output copied to ", dest)
