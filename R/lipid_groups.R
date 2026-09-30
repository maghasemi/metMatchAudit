# Group-name lipid rows, such as
# "sphingomyelin (d18:1/20:2, d18:2/20:1, d16:1/22:2)*".
#
# Metabolon reports these when it could not separate molecular species
# of the same mass. Within one lipid class, the same mass means the same
# chain sum, so every listed member should add up to the same total.
# In the example data all 28 such rows do (notes/findings.md). We use
# that total only when SwissLipids, reached through the row's InChIKey,
# names a Species with the same total. Two independent routes then agree
# on the Species, and the class comes from SwissLipids, not from the
# words in the name.

# Chain sum of one member, e.g. "d18:1/20:2" gives "38:3". Members with
# ether, plasmalogen or hydroxy notation are not handled and give NA,
# because their counting rules differ.
.member_chain_sum <- function(member) {
  if (grepl("O-|P-|\\(OH\\)|;O", member)) return(NA_character_)
  nums <- regmatches(member, gregexpr("\\d+:\\d+", member))[[1]]
  if (length(nums) == 0) return(NA_character_)
  carbons <- sum(as.integer(sub(":.*", "", nums)))
  bonds   <- sum(as.integer(sub(".*:", "", nums)))
  paste0(carbons, ":", bonds)
}

#' Chain sum shared by all members of a group name
#'
#' Reads the comma-separated list of chain notations inside parentheses,
#' adds up the chains of each member, and returns the total if every
#' member gives the same one.
#'
#' @param x A character vector of names.
#'
#' @return A character vector such as `"38:3"`, with `NA` where the name
#'   has no chain list, a member cannot be read, or members disagree.
#'
#' @examples
#' group_chain_sum("sphingomyelin (d18:1/20:2, d18:2/20:1, d16:1/22:2)*")
#' group_chain_sum("diacylglycerol (16:1/18:2 [2], 16:0/18:3 [1])*")
#'
#' @export
group_chain_sum <- function(x) {
  x <- as.character(x)
  vapply(x, function(name) {
    if (is.na(name)) return(NA_character_)
    list_text <- regmatches(name, regexpr("\\([^()]*:[^()]*,[^()]*\\)", name))
    if (length(list_text) == 0) return(NA_character_)
    members <- trimws(strsplit(sub("^\\((.*)\\)$", "\\1", list_text), ",")[[1]])
    sums <- vapply(members, .member_chain_sum, character(1))
    if (anyNA(sums) || length(unique(sums)) != 1) return(NA_character_)
    sums[[1]]
  }, character(1), USE.NAMES = FALSE)
}

# The "carbons:bonds" part of a Species name, e.g. "SM 38:3;O2" gives
# "38:3".
.species_sum <- function(species) {
  out <- regmatches(species, regexpr("\\d+:\\d+", species))
  res <- rep(NA_character_, length(species))
  res[regexpr("\\d+:\\d+", species) > 0] <- out
  res
}
