#' Flag rows with problems in their own evidence
#'
#' Row checks from `notes/design.md`:
#'
#' * C01 `invalid_key`: the InChIKey column holds a value that is not a
#'   standard InChIKey. A missing key is not flagged; many rows simply
#'   have none.
#' * C03 `lipid_source_conflict`: SwissLipids and the name give different
#'   Species for the row, so it got no lipid profile.
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

  has_value <- !is.na(evidence$inchikey) & nzchar(trimws(evidence$inchikey))
  c01 <- has_value & !evidence$ik_valid
  c03 <- !is.na(evidence$lipid_conflict) & evidence$lipid_conflict

  out <- rbind(
    data.frame(row_id = evidence$row_id[c01], check = rep("C01", sum(c01)),
               flag = rep("invalid_key", sum(c01)),
               # sprintf, not paste0: paste0 turns an empty input into one
               # string and would create a phantom row.
               detail = sprintf("InChIKey value: %s", evidence$inchikey[c01]),
               stringsAsFactors = FALSE),
    data.frame(row_id = evidence$row_id[c03], check = rep("C03", sum(c03)),
               flag = rep("lipid_source_conflict", sum(c03)),
               detail = sprintf("SwissLipids Species: %s; the name gives a different Species",
                               evidence$sl_species_abbrev[c03]),
               stringsAsFactors = FALSE)
  )
  rownames(out) <- NULL
  out
}
