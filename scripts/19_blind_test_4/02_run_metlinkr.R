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

# metLinkR's queryRefMet() returns NA for a table with a single column
# (`length(input_df) == 1`), and the run then stops in makeStereochemFlags().
# F2 and F3 have only the name column, so their run copies get one extra,
# empty column that the manifest does not name. The values are unchanged,
# and the copies in data-raw/blind4/ are not touched.
padded <- character(0)
for (f in file.path(run_dir, basename(inputs))) {
  d <- read.csv(f, colClasses = "character", check.names = FALSE, na.strings = character(0))
  if (ncol(d) == 1) {
    d$UNUSED <- ""
    write.csv(d, f, row.names = FALSE)
    padded <- c(padded, basename(f))
  }
}

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

# Pin RaMP DB 3.0.7, the version of every earlier run. metLinkR calls
# RaMP::RaMP() without a version, which asks GitHub for the newest one: on
# 2026-11-08 that request timed out (HTTP 504) and the run stopped, and a
# newer version would otherwise be downloaded and used. With the version
# given and the file already on this machine, GitHub is not asked.
ramp_version <- "3.0.7"
stopifnot(ramp_version %in% RaMP:::listRaMPVersions(local = TRUE))
ramp_orig <- RaMP::RaMP
assignInNamespace("RaMP", function(version = ramp_version, branch = "main") ramp_orig(version, branch), "RaMP")

# metLinkR's pivot_mapping_library() stopped twice inside tidyr's
# names_sep splitting ("SET_VECTOR_ELT() can only be applied to a 'list',
# not a 'promise'"). On the saved input table, the same call ran without
# error in a fresh session, so the failure depends on the state of the long
# run, not on the data. If it happens again, the same pivot is done with
# names_pattern, which avoids that code. On this run's table the two give
# identical results (checked 2026-11-08). The run log says which was used.
pivot_orig <- get("pivot_mapping_library", asNamespace("metLinkR"))
pivot_fallback_used <- FALSE
assignInNamespace("pivot_mapping_library", function(mappings) {
  tryCatch(pivot_orig(mappings), error = function(e) {
    pivot_fallback_used <<- TRUE
    message("pivot_mapping_library() failed (", conditionMessage(e), "); using names_pattern instead")
    mappings |>
      tidyr::pivot_longer(cols = -"Harmonized name", names_to = c(".value", "Origin file"),
                          names_pattern = "^(.*) \\((.*)\\)$", values_drop_na = TRUE) |>
      dplyr::filter(`Input name` != "-") |>
      dplyr::mutate(`Origin file` = gsub("\\)", "", `Origin file`))
  })
}, "metLinkR")

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
  paste("RaMP database version used (pinned):", ramp_version),
  paste("Files given an empty extra column for metLinkR:", paste(padded, collapse = ", ")),
  paste("pivot_mapping_library() replaced by the names_pattern pivot after an error:", pivot_fallback_used),
  "Settings: MetLinkR defaults, n_cores = 1",
  "",
  utils::capture.output(sessionInfo())
), "run_log.txt")
setwd(old)

dir.create(dest, recursive = TRUE)
file.copy(list.files(file.path(run_dir, "metLinkR_output"), full.names = TRUE), dest, overwrite = TRUE)
file.copy(file.path(run_dir, c("manifest.csv", "run_log.txt")), dest, overwrite = TRUE)
message("MetLinkR finished in ", minutes, " minutes. Output copied to ", dest)
