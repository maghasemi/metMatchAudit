# Stereo layers of standard InChIs, and how two of them relate.
#
# Only meaningful for two InChIs with the same main layer (the same
# first 14 InChIKey characters): then the atom numbers in /t and /b refer
# to the same atoms.

#' Read the stereo layers of an InChI
#'
#' Reads the main-layer stereo of a standard InChI:
#'
#' * `/t` tetrahedral centres, for example `t13-,14?,15+` (`?` means
#'   undefined and is dropped)
#' * `/m` `1` means the molecule is the mirror image of what `/t` lists,
#'   so every sign in `/t` flips
#' * `/s` `1` absolute, `2` relative, `3` racemic
#' * `/b` double bonds, for example `b7-6+` (`+` trans, `-` cis)
#'
#' A fixed-H or reconnected layer (`/f`, `/r`) can repeat these; reading
#' stops there.
#'
#' @param inchi One standard InChI string.
#'
#' @return A list with `centres`, a named character vector of `"+"` and
#'   `"-"` (names such as `"t21"` or `"b7-6"`, signs already corrected
#'   for `/m`), and `s`, the `/s` value or `NA`.
#'
#' @export
inchi_stereo <- function(inchi) {
  empty <- list(centres = stats::setNames(character(0), character(0)), s = NA_character_)
  if (length(inchi) != 1 || is.na(inchi) || !startsWith(inchi, "InChI=")) return(empty)
  parts <- strsplit(inchi, "/", fixed = TRUE)[[1]][-1]
  stop_at <- which(substr(parts, 1, 1) %in% c("f", "r"))
  if (length(stop_at) > 0) parts <- parts[seq_len(stop_at[1] - 1)]
  layer <- function(letter) {
    hit <- parts[substr(parts, 1, 1) == letter]
    if (length(hit) == 0) NA_character_ else substring(hit[1], 2)
  }
  read_items <- function(text, prefix, flip) {
    if (is.na(text)) return(character(0))
    items <- strsplit(text, ",", fixed = TRUE)[[1]]
    where <- sub("[-+?]$", "", items)
    sign  <- substring(items, nchar(items))
    if (flip) sign <- chartr("+-", "-+", sign)
    keep <- sign %in% c("+", "-")
    stats::setNames(sign[keep], paste0(prefix, where[keep]))
  }
  m <- layer("m")
  centres <- c(read_items(layer("t"), "t", flip = !is.na(m) && m == "1"),
               read_items(layer("b"), "b", flip = FALSE))
  list(centres = centres, s = layer("s"))
}

#' Compare the stereo of two InChIs with the same connectivity
#'
#' @param inchi_a,inchi_b Standard InChI strings whose InChIKeys share
#'   their first 14 characters.
#'
#' @return One of:
#'   * `"same"`: the same defined centres with the same signs
#'   * `"a_subset_of_b"`, `"b_subset_of_a"`: one defines only some of the
#'     other's centres, and they agree wherever both define one
#'   * `"conflict"`: at least one centre defined in both has different
#'     signs
#'   * `"disjoint"`: each defines centres the other lacks, with no
#'     disagreement
#'   * `"relative"`: either InChI has relative or racemic stereo (`/s2`,
#'     `/s3`), which cannot be compared as absolute
#'   * `"no_stereo"`: neither defines any centre
#'
#' In a test on 317 real structures, removing stereo marks never gave
#' `"conflict"` except for three symmetric molecules such as inositol,
#' and inverting a centre never gave a subset (`notes/decisions.md`,
#' 2026-10-01).
#'
#' @export
compare_stereo <- function(inchi_a, inchi_b) {
  a <- inchi_stereo(inchi_a)
  b <- inchi_stereo(inchi_b)
  if ((!is.na(a$s) && a$s != "1") || (!is.na(b$s) && b$s != "1")) return("relative")
  ca <- a$centres
  cb <- b$centres
  if (length(ca) == 0 && length(cb) == 0) return("no_stereo")
  shared <- intersect(names(ca), names(cb))
  if (any(ca[shared] != cb[shared])) return("conflict")
  if (setequal(names(ca), names(cb))) return("same")
  if (all(names(ca) %in% names(cb))) return("a_subset_of_b")
  if (all(names(cb) %in% names(ca))) return("b_subset_of_a")
  "disjoint"
}
