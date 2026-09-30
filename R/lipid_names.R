# rgoslin, one name at a time. In the feasibility run, one name that
# rgoslin could not parse inside a batch gave an odd error, so every
# call stands alone.

.rgoslin_fields <- c("Grammar", "Level", "Species.Name",
                     "Molecular.Species.Name", "Sn.Position.Name",
                     "Sum.Formula")

# rgoslin levels, coarse to fine, mapped to the four levels this package
# compares at. Class and category are too coarse to use.
.lipid_levels <- c("species", "molecular", "sn_position", "structure")

# Grammars whose output we accept. Checked on every lipid name in the
# MetLinkR example data (notes/decisions.md, 2026-09-30). Shorthand2020,
# Goslin and LipidMaps gave correct results in every case. FattyAcids,
# the grammar for systematic and trivial fatty acid names, did not:
# "3-hydroxyisobutyrate" came out as a wax ester (WE 4:0;O) and
# "Pentadecenoic acid" as FA 15:0 instead of 15:1. It is left out, and
# so is any grammar we have not checked.
.accepted_grammars <- c("Shorthand2020", "Goslin", "LipidMaps")
.rgoslin_level_map <- c(
  SPECIES            = "species",
  MOLECULAR_SPECIES  = "molecular",
  SN_POSITION        = "sn_position",
  STRUCTURE_DEFINED  = "structure",
  FULL_STRUCTURE     = "structure",
  COMPLETE_STRUCTURE = "structure"
)

#' Parse lipid names with rgoslin
#'
#' Each name is parsed on its own. A name counts as parsed only when
#' rgoslin reads it with one of the accepted grammars (Shorthand2020,
#' Goslin, LipidMaps) and reports a level of species or finer. The
#' FattyAcids grammar is not accepted, because on real data it turned
#' small organic acids into wax esters and miscounted double bonds.
#' rgoslin fills some not-applicable fields with the text `"NA"`; these
#' come back as real missing values.
#'
#' @param x A character vector of names.
#'
#' @return A data frame with one row per element of `x`: `parsed`
#'   (logical), `grammar` (the rgoslin grammar, reported even when the
#'   result is not accepted), `level` (one of `"species"`, `"molecular"`,
#'   `"sn_position"`, `"structure"`), `species`, `molecular`,
#'   `sn_position`, `formula`.
#'
#' @examples
#' \donttest{
#' parse_lipid_names(c("SM 36:1;O2", "DG 18:1_18:2", "asparagine"))
#' }
#'
#' @export
parse_lipid_names <- function(x) {
  x <- as.character(x)
  u <- unique(x[!is.na(x) & nzchar(x)])
  parsed <- lapply(u, .parse_one_lipid)
  table <- if (length(parsed) > 0) do.call(rbind, parsed) else .empty_parse(0)
  idx <- match(x, u)
  # Rows for missing or empty names come out of table[NA, ] as all NA.
  out <- table[idx, , drop = FALSE]
  out$parsed[is.na(out$parsed)] <- FALSE
  rownames(out) <- NULL
  out
}

.empty_parse <- function(n) {
  data.frame(parsed = rep(FALSE, n),
             grammar = rep(NA_character_, n),
             level = rep(NA_character_, n),
             species = rep(NA_character_, n),
             molecular = rep(NA_character_, n),
             sn_position = rep(NA_character_, n),
             formula = rep(NA_character_, n),
             stringsAsFactors = FALSE)
}

.fix_na <- function(v) {
  v <- as.character(v)
  v[v %in% c("NA", "")] <- NA_character_
  v
}

.parse_one_lipid <- function(name) {
  res <- NULL
  utils::capture.output(
    res <- tryCatch(suppressMessages(rgoslin::parseLipidNames(name)),
                    error = function(e) NULL)
  )
  out <- .empty_parse(1)
  if (is.null(res) || nrow(res) == 0) return(out)

  field <- function(f) if (f %in% names(res)) .fix_na(res[[f]][1]) else NA_character_
  grammar <- field("Grammar")
  level   <- unname(.rgoslin_level_map[field("Level")])
  out$grammar <- grammar

  if (is.na(grammar) || !grammar %in% .accepted_grammars || is.na(level)) return(out)

  out$parsed      <- TRUE
  out$level       <- level
  out$species     <- field("Species.Name")
  out$molecular   <- field("Molecular.Species.Name")
  out$sn_position <- field("Sn.Position.Name")
  out$formula     <- field("Sum.Formula")
  out
}
