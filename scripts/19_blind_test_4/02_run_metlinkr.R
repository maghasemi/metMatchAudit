# Runs MetLinkR on the three files of the fourth blind test of metMatchAudit
# (notes/decisions.md in that project, 2026-11-08), with the same MetLinkR,
# RaMP and settings as before.
#
# Run from the root of the separate MetLinkR project: Rscript run_metlinkr_blind4.R
# Output goes to metMatchAudit/data-raw/blind4/metlinkr_output/, which git
# ignores there.

suppressPackageStartupMessages({
  library(metLinkR)
  library(RaMP)
})

mma     <- path.expand("~/Projects/metMatchAudit")
inputs  <- file.path(mma, "data-raw/blind4", c("F1.csv", "F2.csv", "F3.csv"))
dest    <- file.path(mma, "data-raw/blind4/metlinkr_output")
run_dir <- "blind4_run"

missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0) stop("Cannot find: ", paste(missing, collapse = ", "))
if (dir.exists(dest)) stop(dest, " exists already. Remove it on purpose to rerun.")

dir.create(run_dir, showWarnings = FALSE)
file.copy(inputs, run_dir, overwrite = TRUE)

# F1 has PubChem and HMDB columns; F2 and F3 have names only.
manifest <- data.frame(
  FileNames       = basename(inputs),
  ShortFileName   = c("F1", "F2", "F3"),
  HMDB            = c("HMDB", NA, NA),
  Metabolite_Name = "BIOCHEMICAL",
  PubChem_CID     = c("PUBCHEM", NA, NA),
  KEGG            = NA,
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
