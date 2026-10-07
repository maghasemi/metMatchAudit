# Runs MetLinkR on the two files of the second blind test of metMatchAudit
# (notes/decisions.md in that project, 2026-10-07), with the same MetLinkR,
# RaMP and settings as the COMETS run.
#
# Run from the root of this project: Rscript run_metlinkr_blind2.R
# Output is copied to metMatchAudit/data-raw/own/metlinkr_output/, which
# git ignores there.

suppressPackageStartupMessages({
  library(metLinkR)
  library(RaMP)
})

src     <- path.expand("~/Projects/metMatchAudit/data-raw/own")
dest    <- file.path(src, "metlinkr_output")
run_dir <- "blind2_run"

files <- c("A_metabolon.csv", "B_ding2021.csv")
missing <- files[!file.exists(file.path(src, files))]
if (length(missing) > 0) stop("Cannot find in ", src, ": ", paste(missing, collapse = ", "))

dir.create(run_dir, showWarnings = FALSE)
file.copy(file.path(src, files), run_dir, overwrite = TRUE)

manifest <- data.frame(
  FileNames       = files,
  ShortFileName   = c("A_metabolon", "B_ding2021"),
  HMDB            = c("HMDB_ID", NA),
  Metabolite_Name = "BIOCHEMICAL",
  PubChem_CID     = "PUBCHEM",
  KEGG            = c("KEGG", NA),
  LIPIDMAPS       = NA,
  chebi           = NA
)
write.csv(manifest, file.path(run_dir, "manifest.csv"), row.names = FALSE, na = "NA")

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
