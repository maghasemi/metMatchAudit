# Expected values are the outputs rgoslin 1.10.0 gave for these names
# in the lipid feasibility run (notes/findings.md).

test_that("shorthand names parse at the right level", {
  out <- parse_lipid_names(c("SM 36:1;O2", "DG 18:1_18:2", "PC 16:0/18:1", "CAR 10:1"))
  expect_true(all(out$parsed))
  expect_identical(out$level, c("species", "molecular", "sn_position", "sn_position"))
  expect_identical(out$species, c("SM 36:1;O2", "DG 36:3", "PC 34:1", "CAR 10:1"))
  expect_identical(out$molecular[2:3], c("DG 18:1_18:2", "PC 16:0_18:1"))
})

test_that("rgoslin's text 'NA' comes back as a real missing value", {
  out <- parse_lipid_names("SM 36:1;O2")
  expect_true(is.na(out$molecular))
})

test_that("Metabolon's long names do not parse", {
  out <- parse_lipid_names("stearoyl sphingomyelin (d18:1/18:0)")
  expect_false(out$parsed)
  expect_true(is.na(out$species))
})

test_that("missing, empty and repeated names are handled", {
  out <- parse_lipid_names(c(NA, "", "SM 36:1;O2", "SM 36:1;O2"))
  expect_equal(nrow(out), 4)
  expect_identical(out$parsed, c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(out$species[3], out$species[4])
})

test_that("an empty input gives an empty result", {
  expect_equal(nrow(parse_lipid_names(character(0))), 0)
})

test_that("the FattyAcids grammar is not accepted", {
  # Real misreadings from the example data: a small organic acid read
  # as a wax ester, and a monounsaturated acid read as saturated.
  out <- parse_lipid_names(c("3-hydroxyisobutyrate", "Pentadecenoic acid"))
  expect_identical(out$parsed, c(FALSE, FALSE))
  expect_identical(out$grammar, c("FattyAcids", "FattyAcids"))
  expect_true(all(is.na(out$species)))
})

test_that("trivial names read by Goslin or LipidMaps grammars are kept", {
  out <- parse_lipid_names(c("Oleic acid", "cholesterol"))
  expect_identical(out$parsed, c(TRUE, TRUE))
  expect_identical(out$species, c("FA 18:1", "ST 27:1;O"))
})

test_that("SwissLipids abbreviations parse into the shared notation", {
  # Checked with real rgoslin on 2026-09-30; all use the LipidMaps grammar.
  out <- parse_lipid_names(c("DG(36:3)", "SM(d36:1)", "HexCer(d34:1)", "DG(18:1_18:2)"))
  expect_identical(out$species, c("DG 36:3", "SM 36:1;O2", "HexCer 34:1;O2", "DG 36:3"))
  expect_identical(out$molecular[4], "DG 18:1_18:2")
  expect_true(all(out$grammar == "LipidMaps"))
})
