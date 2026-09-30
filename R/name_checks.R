# Checks on the text of a metabolite name. Each check is a fixed
# regular expression, so the same name always gives the same answer.
# The patterns were tuned on Metabolon and Broad names from the MetLinkR
# example data. Where a pattern could not tell two meanings apart, it
# was left out on purpose. A missed match makes a later rule fall back
# to "unresolved", which is safe. A false match could make a rule claim
# something the evidence does not support.

# Remove everything inside round or square brackets, including nested
# brackets, so we can see what the name says at its top level.
strip_brackets <- function(x) {
  repeat {
    y <- gsub("\\([^()]*\\)|\\[[^][]*\\]", "", x)
    if (identical(y, x)) return(y)
    x <- y
  }
}

#' Does the name report several compounds at once?
#'
#' Some rows report a group of compounds that the platform could not
#' separate, for example `"(2 or 3)-decenoate"`, `"13-HODE + 9-HODE"`,
#' `"2-Hydroxycaproic and or 2-Hydroxyisocaproic acid"`,
#' `"mannitol/sorbitol"` or `"sphingomyelin (d18:1/20:2, d18:2/20:1)*"`.
#' The InChIKey of such a row describes one member picked from the
#' group, not a confirmed structure.
#'
#' Detected forms:
#'
#' * `" or "` inside round or square brackets: `"(12 or 13)-methylmyristate"`
#' * `" or "` between numbers or locants: `"myo-inositol 1,4 or 1,3-diphosphate"`
#' * `"and or"`, `"and/or"`, `" + "`
#' * two words joined by a slash: `"butyrate/isobutyrate"`
#' * a comma-separated list of lipid chains in parentheses
#' * the word `"Isobar"`
#'
#' A slash between chain notations such as `"18:1/18:2"` is not a group.
#' A plain `" or "` between two full names is not treated as a group,
#' because some files use it to list synonyms of one compound. See
#' [has_alternative_names()].
#'
#' @param x A character vector of names.
#' @return A logical vector the same length as `x`. Missing names give
#'   `FALSE`.
#' @export
is_group_name <- function(x) {
  x <- as.character(x)
  hit <- grepl("\\([^()]* or [^()]*\\)|\\[[^][]* or [^][]*\\]", x, ignore.case = TRUE) |
    grepl("\\d or \\d", x) |
    grepl("and/or|and or| \\+ ", x, ignore.case = TRUE) |
    grepl("[A-Za-z]{3,}/[A-Za-z]{3,}", x) |
    grepl("\\([^()]*:[^()]*,[^()]*\\)", x) |
    grepl("\\bisobar\\b", x, ignore.case = TRUE)
  !is.na(x) & hit
}

#' Does the name list alternative names at its top level?
#'
#' Some files join two names with `" or "`, as in
#' `"Sarcosine or Sarcosine (N-Methylglycine)"` or
#' `"1-palmitoyl-GPC (16:0) or 1-palmitoylglycerophosphocholine (16:0)"`.
#' Often the two sides are synonyms of one compound, but not always:
#' `"branched-chain, straight-chain, or cyclopropyl 10:1 fatty acid"`
#' names real alternatives. The text alone cannot tell these apart, so
#' this check only reports that the name has the form. A row with such a
#' name does not get a confirmed structure from its InChIKey.
#'
#' Only `" or "` outside brackets counts, and not when it sits between
#' numbers (that form is a group, see [is_group_name()]).
#'
#' @inheritParams is_group_name
#' @return A logical vector the same length as `x`.
#' @export
has_alternative_names <- function(x) {
  x <- as.character(x)
  top <- strip_brackets(ifelse(is.na(x), "", x))
  # "and or" and "and/or" belong to is_group_name(), not here.
  top <- gsub("and ?/?or", "", top, ignore.case = TRUE)
  hit <- grepl(" or ", top, ignore.case = TRUE) & !grepl("\\d or \\d", x)
  !is.na(x) & hit
}

#' Is the identification marked as tentative?
#'
#' Metabolon marks compounds that were not confirmed against a
#' standard with a trailing `*`, and with `**` for a lower level of
#' confidence.
#'
#' @inheritParams is_group_name
#' @return A logical vector the same length as `x`.
#' @export
is_tentative_name <- function(x) {
  x <- as.character(x)
  !is.na(x) & grepl("\\*+\\s*$", x)
}

#' Does the name carry an isomer label?
#'
#' Metabolon numbers isomers it measured separately but could not fully
#' assign, as in `"2-methoxyhydroquinone sulfate (1)"` or
#' `"oleoyl-linoleoyl-glycerol (18:1/18:2) [2]"`. Two rows that differ
#' only in this label are different compounds.
#'
#' Only a single digit in round or square brackets at the end of the
#' name counts, after a space and optionally followed by `*`. Without the
#' space rule, peptide names such as `"bradykinin, des-arg(9)"`, where
#' the digit is a residue position, would match.
#'
#' @inheritParams is_group_name
#' @return A logical vector the same length as `x`.
#' @export
has_isomer_label <- function(x) {
  x <- as.character(x)
  !is.na(x) & grepl("\\s(\\(\\d\\)|\\[\\d\\])\\s*\\**\\s*$", x)
}

#' Is the row an internal standard?
#'
#' Internal standards are added to samples for quality control and are
#' not biological measurements. In the MetLinkR example output,
#' `"Inosine-15N4 [iSTD]"` was grouped with inosine.
#'
#' Detected forms: `"[iSTD]"`, `"ISTD"` as a word, `"internal standard"`,
#' and isotope labels such as `13C6`, `15N4`, `18O2`, `2H3`, `-d3` or `d9-`.
#' An isotope label must not follow a letter or digit, so a molecular
#' formula such as `"glucuronide of C10H18O2"` does not match. A `d`
#' followed by a colon, as in the sphingoid base `d18:1`, is not an
#' isotope label either.
#'
#' @inheritParams is_group_name
#' @return A logical vector the same length as `x`.
#' @export
is_internal_standard <- function(x) {
  x <- as.character(x)
  hit <- grepl("\\bistd\\b|internal standard", x, ignore.case = TRUE) |
    grepl("(^|[^A-Za-z0-9])(13C|15N|18O)\\d*", x) |
    grepl("(^|[^A-Za-z0-9])2H\\d+", x) |
    grepl("(^|[-_ (])d\\d{1,2}(?![0-9:])", x, perl = TRUE)
  !is.na(x) & hit
}

#' Does the name state stereochemistry?
#'
#' Used to decide whether a stereo difference between two InChIKeys is
#' backed by the names, or only by keys that may claim more than was
#' measured (see `notes/design.md`, rules R10 and R11).
#'
#' Detected forms: R/S descriptors in parentheses (`"(2S)"`, `"(R)"`,
#' `"(2R,3R)"`), E/Z descriptors in parentheses (`"(E)"`, `"(9Z)"`,
#' `"(2E,4E)"`), `cis-` and `trans-`, `D-` and `L-` before a lowercase
#' letter (`"L-alanine"`), and optical signs `(+)-` and `(-)-`.
#'
#' Not detected, on purpose: `alpha` and `beta`, because in names such
#' as `"alpha-ketoglutarate"` they mark a position, not a configuration.
#' Racemic markers such as `"DL-"` and `"(RS)"` do not count as stated
#' stereo either.
#'
#' @inheritParams is_group_name
#' @return A logical vector the same length as `x`.
#' @export
states_stereo <- function(x) {
  x <- as.character(x)
  hit <- grepl("\\((\\d+[a-z]?)?[RS](,\\s*(\\d+[a-z]?)?[RS])*\\)", x) |
    grepl("\\((\\d+)?[EZ](,\\s*(\\d+)?[EZ])*\\)", x) |
    grepl("\\b(cis|trans)-", x, ignore.case = TRUE) |
    grepl("(^|[^A-Za-z])[DL]-(?=[a-z])", x, perl = TRUE) |
    grepl("\\([+-]\\)-", x)
  !is.na(x) & hit
}
