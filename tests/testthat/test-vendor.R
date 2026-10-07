test_that("the Biocrates table maps lipids and plain acylcarnitines only", {
  tab <- biocrates_p180_table()
  expect_identical(nrow(tab), 188L)
  get <- function(old) tab$mapped[tab$old == old]
  expect_identical(get("PC ae C34:1"), "PC O-34:1")
  expect_identical(get("lysoPC a C28:0"), "PC O-28:0")
  expect_identical(get("SM (OH) C22:1"), "SM 41:1")
  expect_identical(get("C14:2"), "CAR 14:2")
  expect_identical(get("C16:1-OH"), "CAR 16:1;OH")
  expect_true(is.na(get("C0")))
  expect_true(is.na(get("C5-DC")))     # listed with an alternative
  expect_true(is.na(get("C5-M-DC")))
  expect_true(is.na(get("Ala")))
  expect_true(all(!is.na(tab$species[!is.na(tab$mapped)])))
})

test_that("the abbreviation is read from the last pair of parentheses", {
  expect_identical(.vendor_abbreviation(c(
    "Tetradecadienylcarnitine  ( C14:2 )",
    "Glutarylcarnitine (Hydroxyhexanoylcarnitine)  ( C5-DC(C6-OH) )",
    "Phosphatidylcholine acyl-alkyl C34:1  ( PC ae C34:1 )",
    "PC aa C36:4")),
    c("C14:2", "C5-DC(C6-OH)", "PCaeC34:1", "PCaaC36:4"))
})

vendor_rows <- function(names, species = NA_character_) {
  data.frame(row_id = paste0("r", seq_along(names)), original_name = names,
             lipid_species = species, lipid_molecular = NA_character_,
             lipid_level = ifelse(is.na(species), NA_character_, "structure"),
             lipid_source = ifelse(is.na(species), NA_character_, "swisslipids"),
             lipid_conflict = FALSE, stringsAsFactors = FALSE)
}

test_that("vendor rows get a Species-level profile from the table", {
  ev <- vendor_rows(c("Tetradecadienylcarnitine  ( C14:2 )",
                      "Sphingomyeline C18:0  ( SM C18:0 )",
                      "Glutarylcarnitine (Hydroxyhexanoylcarnitine)  ( C5-DC(C6-OH) )",
                      "Octenoylcarnitine  ( C8:1 )"))
  out <- vendor_lipid_profiles(ev, rep("biocrates_p180", 4))
  expect_identical(out$lipid_species[1], "CAR 14:2")
  expect_identical(out$lipid_species[2], "SM 36:1;O2")
  expect_identical(out$lipid_level[1:2], c("species", "species"))
  expect_identical(out$lipid_source[1:2], c("vendor_table", "vendor_table"))
  # An analyte with an alternative, and one the table does not list.
  expect_true(all(is.na(out$lipid_species[3:4])))
})

test_that("rows of other files are left alone, and an unknown vendor stops", {
  ev <- vendor_rows("Tetradecadienylcarnitine  ( C14:2 )")
  expect_identical(vendor_lipid_profiles(ev, NA_character_), ev)
  expect_error(vendor_lipid_profiles(ev, "nightingale"), "Unknown vendor")
})

test_that("a different Species from another source is a conflict", {
  ev <- vendor_rows("Tetradecadienylcarnitine  ( C14:2 )", species = "CAR 14:1")
  out <- vendor_lipid_profiles(ev, "biocrates_p180")
  expect_true(out$lipid_conflict)
  expect_true(is.na(out$lipid_species))
})
