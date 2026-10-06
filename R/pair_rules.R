# The pair rules from notes/design.md, in the order they run. The first
# rule that applies to a pair decides it. Keep this list, the design
# table and the tests in step with each other.

.rule_order <- c("R02", "R01", "R04h", "R04", "R05g", "R05", "R06", "R06b", "R03",
                 "R08", "R08b", "R09", "R09b", "R10", "R11", "R12", "R99")

.relationships <- c("exact", "broader_narrower", "stereo_variant",
                    "distinct", "unresolved")

.needed_evidence <- c("row_id", "ik_valid", "ik_block1", "ik_block2",
                      "ik_stereo_unspecified", "inchikey", "name_group",
                      "name_alternatives", "name_internal_standard",
                      "name_states_stereo", "lipid_species",
                      "lipid_molecular", "lipid_level")

# Glycosylceramides: SwissLipids writes HexCer and Hex2Cer at Species
# level, while names may say GlcCer, GalCer or LacCer. Map the specific
# forms to the generic one so they can be recognised as naming variants.
.generic_headgroup <- function(species) {
  cls  <- sub(" .*$", "", species)
  rest <- sub("^[^ ]+", "", species)
  cls <- sub("^(Glc|Gal)Cer$", "HexCer", cls)
  cls <- sub("^LacCer$", "Hex2Cer", cls)
  paste0(cls, rest)
}

#' Decide how each pair of rows relates
#'
#' Applies the pair rules from `notes/design.md` to every pair of rows
#' that share a group. The first rule that applies decides, and its ID
#' is kept with the result.
#'
#' | rule | relationship | condition |
#' |---|---|---|
#' | R02 | distinct | either row is an internal standard |
#' | R01 | unresolved | either row has no valid InChIKey and no lipid profile |
#' | R04h | unresolved | lipid Species differ only in hexose headgroup naming (HexCer vs GlcCer) |
#' | R04 | distinct | lipid Species differ |
#' | R05g | unresolved | same Species, one row is a group-name row and the other is reported at molecular level or finer; membership is not checked |
#' | R05 | broader_narrower | same lipid Species, reported at different levels |
#' | R06 | distinct | same Species, both molecular or finer, molecular species differ |
#' | R06b | exact | same Species and same level, at species or molecular level, same name at that level |
#' | R03 | exact | identical InChIKeys |
#' | R08 | unresolved | either name reports a group of compounds |
#' | R08b | unresolved | either name lists alternative names |
#' | R09 | broader_narrower | same connectivity, exactly one key without stereo |
#' | R09b | broader_narrower | same connectivity, and the stereo of one row's trusted structure is a subset of the other's (needs [structure_evidence()]) |
#' | R10 | stereo_variant | same connectivity, different stereo, at least one name states stereo |
#' | R11 | unresolved | same connectivity, different stereo, no name states stereo |
#' | R12 | distinct | different connectivity |
#' | R99 | unresolved | none of the above |
#'
#' Rules R04 to R06b apply only when both rows have a lipid profile.
#' Rules R09 to R12 need a valid InChIKey on both rows. R09b also needs
#' trusted structures from [structure_evidence()]; without them it is
#' skipped.
#'
#' @param evidence Rows after [row_evidence()] and [lipid_evidence()].
#' @param pairs Output of [make_pairs()].
#'
#' @return `pairs` with three columns added: `relationship`, `rule_id`
#'   and `evidence`, a short text naming the values the rule looked at.
#'
#' @export
classify_pairs <- function(evidence, pairs) {
  missing <- setdiff(.needed_evidence, names(evidence))
  if (length(missing) > 0) {
    stop("`evidence` is missing ", paste(missing, collapse = ", "),
         ". Run row_evidence() and lipid_evidence() first.", call. = FALSE)
  }
  ia <- match(pairs$row_a, evidence$row_id)
  ib <- match(pairs$row_b, evidence$row_id)
  if (anyNA(ia) || anyNA(ib)) {
    stop("Some row IDs in `pairs` are not in `evidence`.", call. = FALSE)
  }
  A <- evidence[ia, , drop = FALSE]
  B <- evidence[ib, , drop = FALSE]
  n <- nrow(pairs)

  relationship <- rep(NA_character_, n)
  rule_id      <- rep(NA_character_, n)
  note         <- rep(NA_character_, n)

  decide <- function(cond, rel, id, text) {
    cond <- !is.na(cond) & cond & is.na(relationship)
    if (any(cond)) {
      relationship[cond] <<- rel
      rule_id[cond]      <<- id
      note[cond]         <<- if (length(text) == 1) text else text[cond]
    }
  }

  lipid_both <- !is.na(A$lipid_species) & !is.na(B$lipid_species)
  keys_both  <- A$ik_valid & B$ik_valid
  same_sp    <- lipid_both & A$lipid_species == B$lipid_species
  headgroup_only <- lipid_both & !same_sp &
    .generic_headgroup(A$lipid_species) == .generic_headgroup(B$lipid_species)
  sp_text <- paste0("Species ", A$lipid_species, " vs ", B$lipid_species)
  lv_text <- paste0(A$lipid_species, ": level ", A$lipid_level, " vs ", B$lipid_level)
  mol_text <- paste0("molecular ", A$lipid_molecular, " vs ", B$lipid_molecular)
  finer <- c("molecular", "sn_position", "structure")

  decide(A$name_internal_standard | B$name_internal_standard,
         "distinct", "R02", "one row is an internal standard")

  decide((!A$ik_valid & is.na(A$lipid_species)) | (!B$ik_valid & is.na(B$lipid_species)),
         "unresolved", "R01", "a row has neither a valid InChIKey nor a lipid profile")

  decide(headgroup_only, "unresolved", "R04h", paste0(sp_text, " (headgroup naming only)"))

  decide(lipid_both & !same_sp, "distinct", "R04", sp_text)

  decide(same_sp & ((A$lipid_level == "group" & B$lipid_level %in% finer) |
                      (B$lipid_level == "group" & A$lipid_level %in% finer)),
         "unresolved", "R05g", paste0(lv_text, "; group members not compared"))

  decide(same_sp & A$lipid_level != B$lipid_level, "broader_narrower", "R05", lv_text)

  decide(same_sp & A$lipid_level %in% finer & B$lipid_level %in% finer &
           !is.na(A$lipid_molecular) & !is.na(B$lipid_molecular) &
           A$lipid_molecular != B$lipid_molecular,
         "distinct", "R06", mol_text)

  decide(same_sp & A$lipid_level == B$lipid_level &
           (A$lipid_level == "species" |
              (A$lipid_level == "molecular" & !is.na(A$lipid_molecular) &
                 !is.na(B$lipid_molecular) & A$lipid_molecular == B$lipid_molecular)),
         "exact", "R06b", paste0("same ", A$lipid_level, " annotation ", A$lipid_species))

  decide(keys_both & A$inchikey == B$inchikey, "exact", "R03",
         paste0("identical InChIKey ", A$inchikey))

  decide(A$name_group | B$name_group, "unresolved", "R08",
         "a row reports a group of compounds; overlap not checked")

  decide(A$name_alternatives | B$name_alternatives, "unresolved", "R08b",
         "a row lists alternative names")

  same_b1   <- keys_both & A$ik_block1 == B$ik_block1
  diff_b2   <- same_b1 & A$ik_block2 != B$ik_block2
  one_unspec <- xor(A$ik_stereo_unspecified, B$ik_stereo_unspecified)
  both_spec  <- !A$ik_stereo_unspecified & !B$ik_stereo_unspecified
  key_text   <- paste0(A$inchikey, " vs ", B$inchikey)

  decide(diff_b2 & one_unspec, "broader_narrower", "R09",
         paste0("same connectivity, one key without stereo: ", key_text))

  # R09b uses the full InChI, when structure_evidence() has trusted one
  # for both rows. Only a subset result makes a claim; a conflict falls
  # through to R10 and R11 as before.
  if ("st_usable_inchi" %in% names(A)) {
    both_inchi <- !is.na(A$st_usable_inchi) & !is.na(B$st_usable_inchi)
    stereo_rel <- rep(NA_character_, n)
    todo <- which(diff_b2 & both_inchi & is.na(relationship))
    for (k in todo) stereo_rel[k] <- compare_stereo(A$st_usable_inchi[k], B$st_usable_inchi[k])
    decide(stereo_rel %in% c("a_subset_of_b", "b_subset_of_a"), "broader_narrower", "R09b",
           paste0("same connectivity, stereo of one structure is a subset of the other: ", key_text))
  }

  decide(diff_b2 & both_spec & (A$name_states_stereo | B$name_states_stereo),
         "stereo_variant", "R10", paste0("same connectivity, stereo stated in a name: ", key_text))

  decide(diff_b2 & both_spec, "unresolved", "R11",
         paste0("same connectivity, stereo only in the keys: ", key_text))

  decide(keys_both & A$ik_block1 != B$ik_block1, "distinct", "R12",
         paste0("different connectivity: ", key_text))

  decide(rep(TRUE, n), "unresolved", "R99", "no rule applied")

  # Say where a key came from when it was not in the source file, for
  # every rule that looked at keys.
  if ("inchikey_source" %in% names(A)) {
    key_rules <- rule_id %in% c("R03", "R09", "R09b", "R10", "R11", "R12")
    src <- function(s) ifelse(is.na(s) | s == "file", NA_character_, s)
    sa <- src(A$inchikey_source); sb <- src(B$inchikey_source)
    tag <- ifelse(!is.na(sa) & !is.na(sb), paste0(" (keys from ", sa, " and ", sb, ")"),
           ifelse(!is.na(sa), paste0(" (first key from ", sa, ")"),
           ifelse(!is.na(sb), paste0(" (second key from ", sb, ")"), "")))
    note <- ifelse(key_rules, paste0(note, tag), note)
  }

  pairs$relationship <- relationship
  pairs$rule_id      <- rule_id
  pairs$evidence     <- note
  pairs
}
