# Builds evidence for two rows and classifies the pair. Lipid profile
# columns are given directly, so rule tests do not depend on SwissLipids
# or rgoslin.
pair_result <- function(a, b) {
  blank <- list(inchikey = NA_character_, lipid_species = NA_character_,
                lipid_molecular = NA_character_, lipid_level = NA_character_,
                st_usable_inchi = NA_character_)
  a <- utils::modifyList(blank, a)
  b <- utils::modifyList(blank, b)
  rows <- data.frame(
    row_id          = c("a", "b"),
    original_name   = c(a$name, b$name),
    inchikey        = c(a$inchikey, b$inchikey),
    lipid_species   = c(a$lipid_species, b$lipid_species),
    lipid_molecular = c(a$lipid_molecular, b$lipid_molecular),
    lipid_level     = c(a$lipid_level, b$lipid_level),
    st_usable_inchi = c(a$st_usable_inchi, b$st_usable_inchi),
    stringsAsFactors = FALSE
  )
  ev <- row_evidence(rows)
  classify_pairs(ev, data.frame(group_id = "g", row_a = "a", row_b = "b"))
}
