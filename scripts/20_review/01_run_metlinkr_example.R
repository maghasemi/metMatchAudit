# Review, step 1: runs MetLinkR on the example input that ships with it
# (inst/extdata/HarmInputFiles.csv at commit a6aa145), with default
# settings and RaMP DB 3.0.7, and compares the result with the example
# output that ships next to it (inst/extdata/Example_output/).
#
# The shipped example output was not made from this manifest: its log names
# DDLPS.csv and LEOCC_Metabolon_Annotations.csv, which are not shipped, and
# reads the Broad name from column "name", not "Metabolite"
# (notes/decisions.md, 2026-09-23). Only the three files whose names appear
# in both are compared, and differences are expected for Broad.
#
# Like the other *_run_metlinkr.R scripts, this runs in the separate
# MetLinkR project, not in metMatchAudit. From the root of that project:
#   Rscript ~/Projects/metMatchAudit/scripts/20_review/01_run_metlinkr_example.R
# Output: example_run/ in the MetLinkR project. Nothing is written to
# metMatchAudit.

suppressPackageStartupMessages({
  library(metLinkR)
  library(RaMP)
})

mma      <- path.expand("~/Projects/metMatchAudit")
extdata  <- file.path(mma, "data-raw/MetLinkR/inst/extdata")
shipped  <- file.path(extdata, "Example_output")
run_dir  <- "example_run"

if (!dir.exists(extdata)) stop("Cannot find ", extdata, ". Run scripts/01_setup/01_get_data.sh in metMatchAudit first.")
if (dir.exists(run_dir)) stop(run_dir, " exists already. Remove it on purpose to rerun.")

manifest <- read.csv(file.path(extdata, "HarmInputFiles.csv"), colClasses = "character")
inputs   <- file.path(extdata, manifest$FileNames)
missing  <- inputs[!file.exists(inputs)]
if (length(missing) > 0) stop("Cannot find: ", paste(missing, collapse = ", "))

dir.create(run_dir)
invisible(file.copy(c(inputs, file.path(extdata, "HarmInputFiles.csv")), run_dir))

# Same RaMP pin as in scripts/19_blind_test_4/02_run_metlinkr.R.
ramp_version <- "3.0.7"
stopifnot(ramp_version %in% RaMP:::listRaMPVersions(local = TRUE))
ramp_orig <- RaMP::RaMP
assignInNamespace("RaMP", function(version = ramp_version, branch = "main") ramp_orig(version, branch), "RaMP")

# Same pivot fallback as in scripts/19_blind_test_4/02_run_metlinkr.R. The
# first run of this script (2026-11-16) stopped in pivot_mapping_library()
# with the error seen in the fourth blind test ("SET_VECTOR_ELT() can only
# be applied to a 'list', not a 'promise'"). If the original fails, the
# same pivot is done with names_pattern. The input table is then saved as
# pivot_input.rds, so the two pivots can be compared in a fresh session.
pivot_orig <- get("pivot_mapping_library", asNamespace("metLinkR"))
pivot_fallback_used <- FALSE
assignInNamespace("pivot_mapping_library", function(mappings) {
  tryCatch(pivot_orig(mappings), error = function(e) {
    pivot_fallback_used <<- TRUE
    saveRDS(mappings, "pivot_input.rds")
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
result <- harmonizeInputSheets(inputcsv = "HarmInputFiles.csv", n_cores = 1)
minutes <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)

desc <- utils::packageDescription("metLinkR")
writeLines(c(
  paste("Run date:", format(started, "%Y-%m-%d %H:%M %Z")),
  paste("Run time (minutes):", minutes),
  paste("metLinkR version:", desc$Version),
  paste("metLinkR GitHub commit:", if (is.null(desc$RemoteSha)) "unknown" else desc$RemoteSha),
  paste("RaMP package version:", as.character(utils::packageVersion("RaMP"))),
  paste("RaMP database version used (pinned):", ramp_version),
  paste("pivot_mapping_library() replaced by the names_pattern pivot after an error:", pivot_fallback_used),
  "Settings: MetLinkR defaults, n_cores = 1",
  "",
  utils::capture.output(sessionInfo())
), "run_log.txt")
setwd(old)

# Summary of the new output, one line per file.
out_dir <- file.path(run_dir, "metLinkR_output")
read_out <- function(f) read.csv(f, colClasses = "character", check.names = FALSE, na.strings = c("", "NA"))
mapped <- function(x) !is.na(x) & x != "-"

cat("\nMetLinkR finished in", minutes, "minutes. Pivot fallback used:", pivot_fallback_used, "\n\nNew output:\n")
new_files <- sub("\\.csv$", "_metLinkR.csv", manifest$FileNames)
new <- lapply(file.path(out_dir, new_files), read_out)
names(new) <- new_files
print(data.frame(
  file   = new_files,
  rows   = vapply(new, nrow, 1L),
  mapped = vapply(new, function(d) sum(mapped(d$`Standardized name`)), 1L),
  names  = vapply(new, function(d) length(unique(d$`Standardized name`[mapped(d$`Standardized name`)])), 1L),
  row.names = NULL
))

# Names shared across files: MetLinkR's groups.
long <- do.call(rbind, lapply(names(new), function(f) {
  s <- new[[f]]$`Standardized name`
  data.frame(file = f, name = s[mapped(s)])
}))
files_per_name <- tapply(long$file, long$name, function(x) length(unique(x)))
cat("\nStandardized names found in 2 or more files:", sum(files_per_name >= 2), "\n")
print(table(files_per_name, dnn = "number of files"))

# Comparison with the shipped example output, row by row.
cat("\nAgainst the shipped example output (same row order assumed, checked by row count):\n")
both <- intersect(new_files, list.files(shipped))
for (f in both) {
  a <- new[[f]]
  b <- read_out(file.path(shipped, f))
  if (nrow(a) != nrow(b)) {
    cat(f, ": rows differ,", nrow(a), "new against", nrow(b), "shipped\n")
    next
  }
  sa <- a$`Standardized name`
  sb <- b$`Standardized name`
  same <- (is.na(sa) & is.na(sb)) | (!is.na(sa) & !is.na(sb) & sa == sb)
  cat(f, ":", nrow(a), "rows,", sum(same), "with the same standardized name,", sum(!same), "different\n")
  if (any(!same)) {
    show <- head(which(!same), 10)
    print(data.frame(row = show, new = sa[show], shipped = sb[show]))
  }
}
cat("\nLog:", file.path(run_dir, "run_log.txt"), "\n")
