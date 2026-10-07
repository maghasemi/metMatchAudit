# Runs MetLinkR on the three files of the third blind test of metMatchAudit
# (notes/decisions.md in that project, 2026-10-07), with the same MetLinkR,
# RaMP and settings as before.
#
# Run from the root of this project: Rscript run_metlinkr_blind3.R
# Output goes to metMatchAudit/data-raw/own/metlinkr_output_blind3/, which
# git ignores there, so the second test's output stays as it is.

suppressPackageStartupMessages({
  library(metLinkR)
  library(RaMP)
})

mma     <- path.expand("~/Projects/metMatchAudit")
inputs  <- c(file.path(mma, "data-raw/own/A_metabolon.csv"),
             file.path(mma, "data-raw/lai2021/C1_lai_feces.csv"),
             file.path(mma, "data-raw/lai2021/C2_lai_serum.csv"))
dest    <- file.path(mma, "data-raw/own/metlinkr_output_blind3")
run_dir <- "blind3_run"

missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0) stop("Cannot find: ", paste(missing, collapse = ", "))
if (dir.exists(dest)) stop(dest, " exists already. Remove it on purpose to rerun.")

dir.create(run_dir, showWarnings = FALSE)
file.copy(inputs, run_dir, overwrite = TRUE)

manifest <- data.frame(
  FileNames       = basename(inputs),
  ShortFileName   = c("A_metabolon", "C1_lai_feces", "C2_lai_serum"),
  HMDB            = c("HMDB_ID", NA, NA),
  Metabolite_Name = "BIOCHEMICAL",
  PubChem_CID     = "PUBCHEM",
  KEGG            = c("KEGG", NA, NA),
  LIPIDMAPS       = NA,
  chebi           = NA
)
write.csv(manifest, file.path(run_dir, "manifest.csv"), row.names = FALSE, na = "NA")

old <- setwd(run_dir)
started <- Sys.time()
result <- harmonizeInputSheets(inputcsv = "manifest.csv", n_cores = 1)
minutes <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)

desc <- utils::packageDescription("metLinkR")
writeLines(c(
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
), "run_log.txt")
setwd(old)

dir.create(dest, recursive = TRUE)
file.copy(list.files(file.path(run_dir, "metLinkR_output"), full.names = TRUE), dest, overwrite = TRUE)
file.copy(file.path(run_dir, c("manifest.csv", "run_log.txt")), dest, overwrite = TRUE)
message("MetLinkR finished in ", minutes, " minutes. Output copied to ", dest)
