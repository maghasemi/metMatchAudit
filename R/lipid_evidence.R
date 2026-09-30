#' Build a lipid profile for each row
#'
#' Combines two sources of lipid evidence, as described in
#' `notes/design.md`:
#'
#' * SwissLipids, through the row's InChIKey ([swisslipids_lookup()]).
#'   The Species and Molecular subspecies abbreviations found there are
#'   parsed with rgoslin so they use the same notation as other vendors.
#' * rgoslin, through the row's own name ([parse_lipid_names()]).
#'
#' Rows with a group name or alternative names get no profile, because
#' their InChIKey is not a confirmed structure.
#'
#' When both sources give a Species and they differ, the row is flagged
#' in `lipid_conflict` and gets no profile. When both agree, or only the
#' name parses, the profile follows the name, because the name says what
#' was measured and a key can describe more detail than that. When only
#' SwissLipids gives a Species, the level is `"structure"`, since the
#' row was placed through a full structure.
#'
#' @param rows Output of [row_evidence()].
#' @param db An object from [load_swisslipids()].
#'
#' @return `rows` with these columns added: `sl_tier`, `sl_n_hits`,
#'   `sl_fallback`, `sl_species_abbrev`, `lipid_species`,
#'   `lipid_molecular`, `lipid_level`, `lipid_source`, `lipid_conflict`.
#'
#' @export
lipid_evidence <- function(rows, db) {
  needed <- c("original_name", "inchikey", "ik_valid", "name_group",
              "name_alternatives")
  missing <- setdiff(needed, names(rows))
  if (length(missing) > 0) {
    stop("`rows` is missing ", paste(missing, collapse = ", "),
         ". Run row_evidence() first.", call. = FALSE)
  }
  new_cols <- c("sl_tier", "sl_n_hits", "sl_fallback", "sl_species_abbrev",
                "lipid_species", "lipid_molecular", "lipid_level",
                "lipid_source", "lipid_conflict")
  clash <- intersect(new_cols, names(rows))
  if (length(clash) > 0) {
    stop("`rows` already has these columns: ", paste(clash, collapse = ", "),
         ".", call. = FALSE)
  }

  n <- nrow(rows)
  eligible <- !rows$name_group & !rows$name_alternatives

  # SwissLipids, for eligible rows with a valid key.
  sl <- data.frame(sl_tier = rep(NA_integer_, n), sl_n_hits = 0L,
                   sl_species = NA_character_, sl_molecular = NA_character_,
                   sl_n_species = 0L, sl_fallback = FALSE,
                   stringsAsFactors = FALSE)
  use_sl <- eligible & rows$ik_valid
  if (any(use_sl)) {
    sl[use_sl, ] <- swisslipids_lookup(rows$inchikey[use_sl], db)
  }

  # Put SwissLipids abbreviations into rgoslin notation. Only rows with
  # exactly one Species or one Molecular subspecies are parsed; several
  # means the hierarchy does not settle the question.
  one_species <- ifelse(sl$sl_n_species == 1, sl$sl_species, NA_character_)
  one_molecular <- ifelse(!is.na(sl$sl_molecular) & !grepl(";", sl$sl_molecular),
                          sl$sl_molecular, NA_character_)
  sl_species   <- parse_lipid_names(one_species)$species
  sl_molecular <- parse_lipid_names(one_molecular)$molecular

  # rgoslin on the name, for eligible rows.
  nm <- parse_lipid_names(ifelse(eligible, rows$original_name, NA_character_))

  conflict <- !is.na(sl_species) & nm$parsed & !is.na(nm$species) &
    sl_species != nm$species

  from_name <- nm$parsed & !conflict
  from_sl   <- !from_name & !conflict & !is.na(sl_species)

  species <- ifelse(from_name, nm$species, ifelse(from_sl, sl_species, NA_character_))
  level   <- ifelse(from_name, nm$level, ifelse(from_sl, "structure", NA_character_))
  molecular <- ifelse(from_name & nm$level != "species", nm$molecular,
                      ifelse(from_sl, sl_molecular, NA_character_))
  source <- ifelse(from_name, "rgoslin_name",
                   ifelse(from_sl & sl$sl_fallback, "swisslipids_fallback",
                          ifelse(from_sl, "swisslipids", NA_character_)))

  rows$sl_tier           <- sl$sl_tier
  rows$sl_n_hits         <- sl$sl_n_hits
  rows$sl_fallback       <- sl$sl_fallback
  rows$sl_species_abbrev <- sl$sl_species
  rows$lipid_species     <- species
  rows$lipid_molecular   <- molecular
  rows$lipid_level       <- level
  rows$lipid_source      <- source
  rows$lipid_conflict    <- conflict
  rows
}
