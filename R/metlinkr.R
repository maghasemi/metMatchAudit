#' Read MetLinkR output files into one row table
#'
#' MetLinkR writes each input file back out as `<name>_metLinkR.csv`, with
#' one extra column, `Standardized name`, holding the name it assigned.
#' Rows that share a standardized name form a group. This function reads
#' those files and returns the columns the audit needs.
#'
#' MetLinkR does not record which column of each file holds the
#' metabolite name or the InChIKey, so you say it in `files`.
#'
#' @param files A data frame with one line per MetLinkR output file and
#'   these columns:
#'   * `path`: path to the `_metLinkR.csv` file
#'   * `name_col`: column with the metabolite name
#'   * `inchikey_col`: column with the InChIKey, or `NA` if the file has
#'     none
#'   * `label` (optional): short name used in row IDs; defaults to the
#'     file name without `_metLinkR.csv`
#'   * `smiles_col`, `inchi_col` (optional): columns with a SMILES or an
#'     InChI, used by [structure_evidence()]. `NA` where a file has none.
#'   * `pubchem_col`, `hmdb_col` (optional): columns with PubChem CIDs or
#'     HMDB IDs, used by [resolve_identifiers()].
#' @param encoding Passed to [data.table::fread()]: `"UTF-8"` or
#'   `"Latin-1"`. The MetLinkR example files are Latin-1.
#'
#' @return A data frame with one line per input row: `row_id`
#'   (`label:row number`), `file`, `row_in_file`, `original_name`,
#'   `inchikey`, `smiles`, `inchi`, `pubchem`, `hmdb` and `group_id` (the standardized name;
#'   missing when MetLinkR assigned none).
#'
#' @export
read_metlinkr <- function(files, encoding = "UTF-8") {
  if (!is.data.frame(files)) stop("`files` must be a data frame.", call. = FALSE)
  missing <- setdiff(c("path", "name_col", "inchikey_col"), names(files))
  if (length(missing) > 0) {
    stop("`files` needs the columns ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  for (opt in c("smiles_col", "inchi_col", "pubchem_col", "hmdb_col")) {
    if (!opt %in% names(files)) files[[opt]] <- NA_character_
  }
  if (!"label" %in% names(files)) {
    files$label <- sub("_metLinkR\\.csv$", "", basename(files$path))
  }
  if (anyDuplicated(files$label)) {
    stop("File labels must be unique: ", paste(files$label, collapse = ", "), call. = FALSE)
  }

  pieces <- lapply(seq_len(nrow(files)), function(i) {
    f <- files[i, ]
    if (!file.exists(f$path)) stop("Cannot find ", f$path, ".", call. = FALSE)
    d <- data.table::fread(f$path, colClasses = "character", encoding = encoding,
                           na.strings = c("", "NA"), data.table = FALSE,
                           showProgress = FALSE)
    need <- c(f$name_col, "Standardized name",
              if (!is.na(f$inchikey_col)) f$inchikey_col,
              if (!is.na(f$smiles_col)) f$smiles_col,
              if (!is.na(f$inchi_col)) f$inchi_col,
              if (!is.na(f$pubchem_col)) f$pubchem_col,
              if (!is.na(f$hmdb_col)) f$hmdb_col)
    absent <- setdiff(need, names(d))
    if (length(absent) > 0) {
      stop(f$path, " has no column ", paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    group <- trimws(d[["Standardized name"]])
    group[group %in% c("", "-")] <- NA_character_
    data.frame(
      row_id        = paste0(f$label, ":", seq_len(nrow(d))),
      file          = f$label,
      row_in_file   = seq_len(nrow(d)),
      original_name = d[[f$name_col]],
      inchikey      = if (is.na(f$inchikey_col)) NA_character_ else trimws(d[[f$inchikey_col]]),
      smiles        = if (is.na(f$smiles_col)) NA_character_ else trimws(d[[f$smiles_col]]),
      inchi         = if (is.na(f$inchi_col)) NA_character_ else trimws(d[[f$inchi_col]]),
      pubchem       = if (is.na(f$pubchem_col)) NA_character_ else trimws(d[[f$pubchem_col]]),
      hmdb          = if (is.na(f$hmdb_col)) NA_character_ else trimws(d[[f$hmdb_col]]),
      group_id      = group,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

#' Audit MetLinkR groupings
#'
#' Runs the whole audit: reads the MetLinkR output files, collects the
#' evidence for every row, classifies every pair of rows that MetLinkR
#' put in the same group, flags rows with problems, and summarises each
#' group.
#'
#' @inheritParams read_metlinkr
#' @param swisslipids Either a path to SwissLipids `lipids.tsv` or an
#'   object from [load_swisslipids()]. Loading the file takes about half
#'   a minute, so pass a loaded object when you run several audits.
#' @param expected_md5 Optional md5 of the SwissLipids file, checked when
#'   `swisslipids` is a path.
#' @param use_rdkit Use RDKit for the structure checks (R09b, C04) when
#'   it is available. See [structure_evidence()].
#' @param resolve_ids Fill in missing InChIKeys from PubChem CIDs or HMDB
#'   IDs. See [resolve_identifiers()].
#' @param cache_file Where to keep PubChem answers between runs.
#' @param refmet Optional. A table from [load_refmet()] or the path to the
#'   RefMet download. When given, each row is checked against the RefMet
#'   entry of its group with [refmet_check()], which adds row flags C06
#'   and C07. Pair results do not change.
#' @param refmet_md5 Expected md5 of the RefMet file, when `refmet` is a
#'   path.
#'
#' @return An object of class `metmatch_audit`, a list with:
#'   * `pairs`: every pair with `relationship`, `rule_id`, `evidence`
#'   * `groups`: one line per group, see [summarise_groups()]
#'   * `row_flags`: see [row_flags()]
#'   * `lipid_matches`: lipid pairs across files, see
#'     [propose_lipid_matches()]
#'   * `rows`: every row with its evidence
#'   * `info`: package version, date, and the md5 of the SwissLipids file
#'
#' @export
audit_metlinkr <- function(files, swisslipids, encoding = "UTF-8",
                           expected_md5 = NULL, use_rdkit = TRUE,
                           resolve_ids = TRUE, cache_file = NULL,
                           refmet = NULL, refmet_md5 = NULL) {
  db <- if (inherits(swisslipids, "swisslipids_db")) {
    swisslipids
  } else {
    load_swisslipids(swisslipids, expected_md5 = expected_md5)
  }

  rows     <- read_metlinkr(files, encoding = encoding)
  rows     <- if (resolve_ids) {
    resolve_identifiers(rows, cache_file = cache_file)
  } else {
    rows$inchikey_file <- rows$inchikey
    rows$inchikey_source <- ifelse(is_valid_inchikey(rows$inchikey), "file", NA_character_)
    rows
  }
  evidence <- lipid_evidence(row_evidence(rows), db)
  evidence <- structure_evidence(evidence, use_rdkit = use_rdkit)
  pairs    <- classify_pairs(evidence, make_pairs(evidence))
  if (!is.null(refmet)) {
    if (!inherits(refmet, "refmet_db")) refmet <- load_refmet(refmet, expected_md5 = refmet_md5)
    evidence <- refmet_check(evidence, refmet, cache_file = cache_file)
  }

  out <- list(
    pairs     = pairs,
    lipid_matches = propose_lipid_matches(evidence, group_col = "group_id"),
    groups    = summarise_groups(pairs),
    row_flags = row_flags(evidence),
    rows      = evidence,
    info      = list(
      package_version = as.character(utils::packageVersion("metMatchAudit")),
      run_date        = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      swisslipids_md5 = db$md5,
      refmet_md5      = if (is.null(refmet)) NA_character_ else attr(refmet, "md5"),
      structure_checks = use_rdkit && rdkit_available(),
      keys_by_source = table(rows$inchikey_source, useNA = "ifany"),
      files           = files
    )
  )
  class(out) <- "metmatch_audit"
  out
}

#' @export
print.metmatch_audit <- function(x, ...) {
  cat("<metmatch_audit>\n")
  cat(" rows:  ", nrow(x$rows), " in ", length(unique(x$rows$file)), " files\n", sep = "")
  cat(" groups:", nrow(x$groups), " with two or more rows\n")
  cat(" pairs: ", nrow(x$pairs), "\n", sep = "")
  if (nrow(x$pairs) > 0) {
    tab <- table(factor(x$pairs$relationship, levels = .severity))
    for (r in names(tab)) cat("   ", format(r, width = 17), tab[[r]], "\n")
  }
  cat(" row flags:", nrow(x$row_flags), "\n")
  if (!is.null(x$lipid_matches)) {
    lm <- x$lipid_matches
    cat(" lipid matches across files:", nrow(lm),
        "(", sum(!lm$already_grouped), "not grouped by MetLinkR )\n")
  }
  invisible(x)
}

#' Write an audit to CSV files
#'
#' Writes `pairs.csv`, `groups.csv`, `row_flags.csv`, `rows.csv`,
#' `lipid_matches.csv` and a short `run_info.txt` into `dir`.
#'
#' @param audit An object from [audit_metlinkr()].
#' @param dir Output folder. Created if it does not exist.
#'
#' @return The paths of the written files, invisibly.
#'
#' @export
write_audit <- function(audit, dir) {
  if (!inherits(audit, "metmatch_audit")) {
    stop("`audit` must come from audit_metlinkr().", call. = FALSE)
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- file.path(dir, c("pairs.csv", "groups.csv", "row_flags.csv", "rows.csv",
                            "lipid_matches.csv", "run_info.txt"))
  utils::write.csv(audit$pairs,     paths[1], row.names = FALSE)
  utils::write.csv(audit$groups,    paths[2], row.names = FALSE)
  utils::write.csv(audit$row_flags, paths[3], row.names = FALSE)
  utils::write.csv(audit$rows,      paths[4], row.names = FALSE)
  utils::write.csv(audit$lipid_matches, paths[5], row.names = FALSE)
  info <- audit$info
  writeLines(c(
    paste("metMatchAudit version:", info$package_version),
    paste("Run date:", info$run_date),
    paste("SwissLipids md5:", info$swisslipids_md5),
    paste("RefMet md5:", if (is.null(info$refmet_md5)) NA else info$refmet_md5),
    paste("Structure checks (RDKit):", isTRUE(info$structure_checks)),
    paste("InChIKeys by source:", paste(names(info$keys_by_source), info$keys_by_source,
                                          sep = " ", collapse = ", ")),
    "Input files:",
    paste0("  ", info$files$path, " (name: ", info$files$name_col,
           ", InChIKey: ", info$files$inchikey_col, ")")
  ), paths[6])
  invisible(paths)
}
