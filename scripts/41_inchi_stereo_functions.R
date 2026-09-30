# Stereo layers of a standard InChI, and how two of them relate.
#
# Only meaningful for two InChIs with the same main layer (same first
# 14 characters of the InChIKey): then atom numbers in /t and /b refer
# to the same atoms.
#
#   /t  tetrahedral centres, e.g. t13-,14?,15+   ("?" = undefined)
#   /m  0 or 1; 1 means the actual molecule is the mirror image of what
#       /t lists, so every +/- in /t flips
#   /s  1 absolute, 2 relative, 3 racemic
#   /b  double bonds, e.g. b7-6+,10-9-  ("+" trans, "-" cis, "?" undefined)
#
# Only the main layers are read. A fixed-H or reconnected layer (/f, /r)
# can repeat /t or /b; everything after it is ignored.

inchi_stereo <- function(inchi) {
  empty <- list(centres = character(0), s = NA_character_)
  if (is.na(inchi) || !startsWith(inchi, "InChI=")) return(empty)
  parts <- strsplit(inchi, "/", fixed = TRUE)[[1]][-1]
  stop_at <- which(substr(parts, 1, 1) %in% c("f", "r"))
  if (length(stop_at) > 0) parts <- parts[seq_len(stop_at[1] - 1)]
  layer <- function(letter) {
    hit <- parts[substr(parts, 1, 1) == letter]
    if (length(hit) == 0) NA_character_ else substring(hit[1], 2)
  }
  t_layer <- layer("t"); m_layer <- layer("m"); s_layer <- layer("s"); b_layer <- layer("b")
  out <- character(0)
  if (!is.na(t_layer)) {
    items <- strsplit(t_layer, ",", fixed = TRUE)[[1]]
    atom <- sub("[-+?]$", "", items)
    sign <- substring(items, nchar(items))
    if (!is.na(m_layer) && m_layer == "1") {
      sign <- chartr("+-", "-+", sign)
    }
    keep <- sign %in% c("+", "-")
    out <- c(out, stats::setNames(sign[keep], paste0("t", atom[keep])))
  }
  if (!is.na(b_layer)) {
    items <- strsplit(b_layer, ",", fixed = TRUE)[[1]]
    bond <- sub("[-+?]$", "", items)
    sign <- substring(items, nchar(items))
    keep <- sign %in% c("+", "-")
    out <- c(out, stats::setNames(sign[keep], paste0("b", bond[keep])))
  }
  list(centres = out, s = s_layer)
}

# "same", "a_subset_of_b", "b_subset_of_a", "conflict", "relative" (a /s
# layer other than 1 on either side), or "no_stereo" (neither defines any).
compare_stereo <- function(inchi_a, inchi_b) {
  a <- inchi_stereo(inchi_a); b <- inchi_stereo(inchi_b)
  if ((!is.na(a$s) && a$s != "1") || (!is.na(b$s) && b$s != "1")) return("relative")
  ca <- a$centres; cb <- b$centres
  if (length(ca) == 0 && length(cb) == 0) return("no_stereo")
  shared <- intersect(names(ca), names(cb))
  if (any(ca[shared] != cb[shared])) return("conflict")
  if (setequal(names(ca), names(cb))) return("same")
  if (all(names(ca) %in% names(cb))) return("a_subset_of_b")
  if (all(names(cb) %in% names(ca))) return("b_subset_of_a")
  "conflict_sets"  # each defines centres the other lacks, and they agree where both do
}
