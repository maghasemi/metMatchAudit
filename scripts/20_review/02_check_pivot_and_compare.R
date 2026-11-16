# Review, step 1, after 01_run_metlinkr_example.R. Two checks:
#
# 1. The pivot fallback. In the run of 2026-11-16, pivot_mapping_library()
#    failed again and the names_pattern pivot was used. Here both are run on
#    the saved input table in a fresh session, and the results are compared.
# 2. The new output against the shipped example output, matched on the
#    input columns instead of on row position. Row position does not work:
#    the Broad files have the same row count but a different row order.
#
# Run from the root of the MetLinkR project, after 01_run_metlinkr_example.R:
#   Rscript ~/Projects/metMatchAudit/scripts/20_review/02_check_pivot_and_compare.R

suppressPackageStartupMessages(library(metLinkR))

run_dir <- "example_run"
shipped <- path.expand("~/Projects/metMatchAudit/data-raw/MetLinkR/inst/extdata/Example_output")

# 1. Pivot
f <- file.path(run_dir, "pivot_input.rds")
if (file.exists(f)) {
  m <- readRDS(f)
  a <- tryCatch(metLinkR:::pivot_mapping_library(m), error = function(e) e)
  b <- m |>
    tidyr::pivot_longer(cols = -"Harmonized name", names_to = c(".value", "Origin file"),
                        names_pattern = "^(.*) \\((.*)\\)$", values_drop_na = TRUE) |>
    dplyr::filter(`Input name` != "-") |>
    dplyr::mutate(`Origin file` = gsub("\\)", "", `Origin file`))
  if (inherits(a, "error")) {
    cat("Original pivot fails in a fresh session too:", conditionMessage(a), "\n")
  } else {
    cat("Original pivot in a fresh session:", nrow(a), "rows; fallback:", nrow(b), "rows\n")
    cat("Identical:", isTRUE(all.equal(as.data.frame(a), as.data.frame(b), check.attributes = FALSE)), "\n")
  }
} else {
  cat("No pivot_input.rds: the fallback was not used.\n")
}

# 2. Keyed comparison
read_out <- function(p) read.csv(p, colClasses = "character", check.names = FALSE, na.strings = c("", "NA"))
for (fn in intersect(list.files(file.path(run_dir, "metLinkR_output"), pattern = "_metLinkR\\.csv$"), list.files(shipped))) {
  a <- read_out(file.path(run_dir, "metLinkR_output", fn))
  b <- read_out(file.path(shipped, fn))
  cols <- setdiff(intersect(names(a), names(b)), "Standardized name")
  key  <- function(d) do.call(paste, c(d[cols], sep = "\r"))
  a$key <- key(a)
  b$key <- key(b)
  dup <- sum(duplicated(a$key)) + sum(duplicated(b$key))
  j <- merge(a[!duplicated(a$key), c("key", "Standardized name")],
             b[!duplicated(b$key), c("key", "Standardized name")],
             by = "key", suffixes = c(".new", ".shipped"))
  sa <- j$`Standardized name.new`
  sb <- j$`Standardized name.shipped`
  same <- (is.na(sa) & is.na(sb)) | (!is.na(sa) & !is.na(sb) & sa == sb)
  cat("\n", fn, ": ", nrow(a), " new rows, ", nrow(b), " shipped, ", nrow(j), " matched on ",
      length(cols), " input columns (", dup, " duplicate keys left out); ",
      sum(same), " same standardized name, ", sum(!same), " different\n", sep = "")
  if (any(!same)) {
    d <- data.frame(new = sa[!same], shipped = sb[!same])
    cat("  new empty, shipped has a name:", sum(is.na(d$new) & !is.na(d$shipped)),
        "| new has a name, shipped empty:", sum(!is.na(d$new) & is.na(d$shipped)),
        "| both named, different:", sum(!is.na(d$new) & !is.na(d$shipped)), "\n")
    print(head(d[!is.na(d$new) & !is.na(d$shipped), ], 10), row.names = FALSE)
  }
}
