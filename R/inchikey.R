#' Check whether values are standard InChIKeys
#'
#' A standard InChIKey has 27 characters: 14 capital letters, a hyphen,
#' 10 capital letters, a hyphen, and one capital letter. Anything else
#' returns `FALSE`, including missing values, lowercase text, and keys
#' that still carry an `InChIKey=` prefix. Values are never corrected.
#'
#' @param x A character vector.
#'
#' @return A logical vector the same length as `x`, with no missing
#'   values.
#'
#' @examples
#' is_valid_inchikey(c("BPYKTIZUTYGOLE-KDUUSRDASA-N", "inchikey", NA))
#'
#' @export
is_valid_inchikey <- function(x) {
  x <- as.character(x)
  !is.na(x) & grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", x)
}

#' Split InChIKeys into their parts
#'
#' A standard InChIKey such as `BPYKTIZUTYGOLE-KDUUSRDASA-N` has three
#' parts:
#'
#' * block1, the first 14 letters, encodes the connectivity skeleton.
#'   Two keys with the same block1 describe molecules whose atoms are
#'   connected the same way.
#' * block2, the next 10 letters, encodes stereochemistry, isotopes and
#'   version flags.
#' * the last letter encodes the protonation state.
#'
#' `stereo_unspecified` is `TRUE` when block2 is `UHFFFAOYSA`, which is
#' what a standard InChIKey carries when the structure has no stereo
#' layer and no isotope layer.
#'
#' Invalid keys (see [is_valid_inchikey()]) give `NA` in every column.
#'
#' @param x A character vector of InChIKeys.
#'
#' @return A data frame with one row per element of `x` and columns
#'   `block1`, `block2`, `protonation` and `stereo_unspecified`.
#'
#' @examples
#' inchikey_parts(c("BPYKTIZUTYGOLE-KDUUSRDASA-N",
#'                  "OQWOHRPOYAVIOK-UHFFFAOYSA-N",
#'                  "inchikey"))
#'
#' @export
inchikey_parts <- function(x) {
  x <- as.character(x)
  ok <- is_valid_inchikey(x)
  block1 <- ifelse(ok, substr(x, 1, 14), NA_character_)
  block2 <- ifelse(ok, substr(x, 16, 25), NA_character_)
  data.frame(
    block1             = block1,
    block2             = block2,
    protonation        = ifelse(ok, substr(x, 27, 27), NA_character_),
    stereo_unspecified = ifelse(ok, block2 == "UHFFFAOYSA", NA),
    stringsAsFactors   = FALSE
  )
}
