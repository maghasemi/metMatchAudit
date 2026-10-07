#' Flag rows with problems in their own evidence
#'
#' Row checks from `notes/design.md`:
#'
#' * C01 `invalid_key`: the InChIKey column holds a value that is not a
#'   standard InChIKey. A missing key is not flagged; many rows simply
#'   have none.
#' * C03 `lipid_source_conflict`: SwissLipids and the name give different
#'   Species for the row, so it got no lipid profile.
#' * C04 `structure_key_connectivity`: the row's own structure (InChI or
#'   SMILES) and its stored InChIKey describe molecules with different
#'   connectivity. Needs [structure_evidence()].
#' * C04 `structure_key_stereo`: same connectivity, but the stereo or
#'   isotope block differs.
#' * C06 `refmet_assignment_conflict`: the row's own evidence says it is
#'   not the RefMet entry its group stands for. Needs [refmet_check()].
#' * C07 `refmet_identifier_conflict`: the row and its RefMet entry have
#'   the same name but different structures, so one of the identifiers is
#'   wrong or the two describe different forms. Needs [refmet_check()].
#'
#' @param evidence Rows after [row_evidence()] and [lipid_evidence()].
#'
#' @return A data frame with one line per flag: `row_id`, `check`,
#'   `flag` and `detail`. Rows without problems do not appear.
#'
#' @export
row_flags <- function(evidence) {
  needed <- c("row_id", "inchikey", "ik_valid", "lipid_conflict", "sl_species_abbrev")
  missing <- setdiff(needed, names(evidence))
  if (length(missing) > 0) {
    stop("`evidence` is missing ", paste(missing, collapse = ", "),
         ". Run row_evidence() and lipid_evidence() first.", call. = FALSE)
  }

  # C01 is about the value in the source file. After identifier
  # resolution that value is kept in inchikey_file.
  file_key <- if ("inchikey_file" %in% names(evidence)) evidence$inchikey_file else evidence$inchikey
  has_value <- !is.na(file_key) & nzchar(trimws(file_key))
  c01 <- has_value & !is_valid_inchikey(file_key)
  c03 <- !is.na(evidence$lipid_conflict) & evidence$lipid_conflict

  c04_b1 <- c04_st <- rep(FALSE, nrow(evidence))
  if (all(c("st_key_agrees", "st_block1_differs", "st_key") %in% names(evidence))) {
    c04_b1 <- !is.na(evidence$st_block1_differs) & evidence$st_block1_differs
    c04_st <- !is.na(evidence$st_key_agrees) & !evidence$st_key_agrees & !c04_b1
  }
  c06 <- c07 <- rep(FALSE, nrow(evidence))
  if (all(c("rm_check", "rm_entry", "rm_detail") %in% names(evidence))) {
    c06 <- !is.na(evidence$rm_check) & evidence$rm_check == "conflict"
    c07 <- !is.na(evidence$rm_check) & evidence$rm_check == "identifier_conflict"
  }
  rm_detail <- function(sel) sprintf("RefMet entry %s; %s", evidence$rm_entry[sel], evidence$rm_detail[sel])

  c04_detail <- function(sel) sprintf("stored key %s, key from the row's structure %s",
                                      evidence$inchikey[sel], evidence$st_key[sel])

  out <- rbind(
    data.frame(row_id = evidence$row_id[c01], check = rep("C01", sum(c01)),
               flag = rep("invalid_key", sum(c01)),
               # sprintf, not paste0: paste0 turns an empty input into one
               # string and would create a phantom row.
               detail = sprintf("InChIKey value: %s", file_key[c01]),
               stringsAsFactors = FALSE),
    data.frame(row_id = evidence$row_id[c03], check = rep("C03", sum(c03)),
               flag = rep("lipid_source_conflict", sum(c03)),
               detail = sprintf("SwissLipids Species: %s; the name gives a different Species",
                               evidence$sl_species_abbrev[c03]),
               stringsAsFactors = FALSE),
    data.frame(row_id = evidence$row_id[c04_b1], check = rep("C04", sum(c04_b1)),
               flag = rep("structure_key_connectivity", sum(c04_b1)),
               detail = c04_detail(c04_b1), stringsAsFactors = FALSE),
    data.frame(row_id = evidence$row_id[c04_st], check = rep("C04", sum(c04_st)),
               flag = rep("structure_key_stereo", sum(c04_st)),
               detail = c04_detail(c04_st), stringsAsFactors = FALSE),
    data.frame(row_id = evidence$row_id[c06], check = rep("C06", sum(c06)),
               flag = rep("refmet_assignment_conflict", sum(c06)),
               detail = rm_detail(c06), stringsAsFactors = FALSE),
    data.frame(row_id = evidence$row_id[c07], check = rep("C07", sum(c07)),
               flag = rep("refmet_identifier_conflict", sum(c07)),
               detail = rm_detail(c07), stringsAsFactors = FALSE)
  )
  rownames(out) <- NULL
  out
}
