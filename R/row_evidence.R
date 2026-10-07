#' Collect evidence for each row
#'
#' Adds the InChIKey parts and the name checks to a table of rows. This
#' is the non-lipid part of the evidence described in `notes/design.md`.
#' Nothing in the input is changed; new columns are added on the right.
#'
#' @param rows A data frame with a character column `original_name`.
#'   A column `inchikey` is used if present; without it every key is
#'   treated as missing.
#'
#' @return `rows` with these columns added:
#'   * `ik_valid`, `ik_block1`, `ik_block2`, `ik_protonation`,
#'     `ik_stereo_unspecified`: see [is_valid_inchikey()] and
#'     [inchikey_parts()]
#'   * `name_group`: [is_group_name()]
#'   * `name_alternatives`: [has_alternative_names()]
#'   * `name_tentative`: [is_tentative_name()]
#'   * `name_isomer_label`: [has_isomer_label()]
#'   * `name_isomer_marker`: [has_isomer_marker()]
#'   * `name_internal_standard`: [is_internal_standard()]
#'   * `name_states_stereo`: [states_stereo()]
#'   * `name_normalized`: [normalize_name()]
#'   * `name_stereo_terms`: [stereo_terms()]
#'   * `name_omega`: [omega_positions()]
#'
#' @examples
#' rows <- data.frame(
#'   original_name = c("bilirubin (E,E)*", "Inosine-15N4 [iSTD]"),
#'   inchikey      = c("BPYKTIZUTYGOLE-KDUUSRDASA-N", "inchikey")
#' )
#' row_evidence(rows)
#'
#' @export
row_evidence <- function(rows) {
  if (!is.data.frame(rows)) {
    stop("`rows` must be a data frame.", call. = FALSE)
  }
  if (!"original_name" %in% names(rows)) {
    stop("`rows` needs a column called `original_name`.", call. = FALSE)
  }

  new_cols <- c("ik_valid", "ik_block1", "ik_block2", "ik_protonation",
                "ik_stereo_unspecified", "name_group", "name_alternatives",
                "name_tentative", "name_isomer_label", "name_isomer_marker",
                "name_internal_standard", "name_states_stereo",
                "name_normalized", "name_stereo_terms", "name_omega")
  clash <- intersect(new_cols, names(rows))
  if (length(clash) > 0) {
    stop("`rows` already has these columns: ", paste(clash, collapse = ", "),
         ". Remove them first so nothing is overwritten by accident.",
         call. = FALSE)
  }

  keys  <- if ("inchikey" %in% names(rows)) as.character(rows$inchikey) else rep(NA_character_, nrow(rows))
  names <- as.character(rows$original_name)
  parts <- inchikey_parts(keys)

  rows$ik_valid               <- is_valid_inchikey(keys)
  rows$ik_block1              <- parts$block1
  rows$ik_block2              <- parts$block2
  rows$ik_protonation         <- parts$protonation
  rows$ik_stereo_unspecified  <- parts$stereo_unspecified
  rows$name_group             <- is_group_name(names)
  rows$name_alternatives      <- has_alternative_names(names)
  rows$name_tentative         <- is_tentative_name(names)
  rows$name_isomer_label      <- has_isomer_label(names)
  rows$name_isomer_marker     <- has_isomer_marker(names)
  rows$name_internal_standard <- is_internal_standard(names)
  rows$name_states_stereo     <- states_stereo(names)
  rows$name_normalized        <- normalize_name(names)
  rows$name_stereo_terms      <- stereo_terms(names)
  rows$name_omega             <- omega_positions(names)
  rows
}
