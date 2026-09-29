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
