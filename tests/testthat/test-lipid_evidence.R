make_rows <- function() {
  data.frame(
    row_id = c("Broad:SM", "LEOCC:SM", "LEOCC:DG1", "LEOCC:MG",
               "LEOCC:SMgroup", "LEOCC:CAR", "made-up:conflict", "COMETS:sarc"),
    original_name = c(
      "SM 36:1;O2",
      "stearoyl sphingomyelin (d18:1/18:0)",
      "oleoyl-linoleoyl-glycerol (18:1/18:2) [1]",
      "1-oleoylglycerol (18:1)",
      "sphingomyelin (d18:1/20:2, d18:2/20:1)*",
      "cis-4-decenoylcarnitine (C10:1)",
      "DG 18:1_18:2",
      "Sarcosine or Sarcosine (N-Methylglycine)"
    ),
    inchikey = c(
      "LKQLRGMMMAHREN-YJFXYUILSA-N",
      "LKQLRGMMMAHREN-YJFXYUILSA-N",
      "BLZVZPYMHLXLHG-RQOIEFAZSA-N",
      "RZRNAYUHWVFMIP-KTKRTIGZSA-N",
      "BLZVZPYMHLXLHG-RQOIEFAZSA-N",
      "OQWOHRPOYAVIOK-KTKRTIGZSA-N",
      "LKQLRGMMMAHREN-YJFXYUILSA-N",
      NA
    ),
    stringsAsFactors = FALSE
  )
}

evidence <- function() {
  db <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  lipid_evidence(row_evidence(make_rows()), db)
}

test_that("Broad and Metabolon names for one SM reach the same Species", {
  out <- evidence()
  expect_identical(out$lipid_species[1:2], c("SM 36:1;O2", "SM 36:1;O2"))
  # The Broad name reports a Species; the Metabolon row was placed
  # through its structure.
  expect_identical(out$lipid_level[1:2], c("species", "structure"))
  expect_identical(out$lipid_source[1:2], c("rgoslin_name", "swisslipids"))
})

test_that("a structure row gets Species and Molecular species from SwissLipids", {
  out <- evidence()
  expect_identical(out$lipid_species[3], "DG 36:3")
  expect_identical(out$lipid_molecular[3], "DG 18:1_18:2")
})

test_that("the connectivity fallback is marked as such", {
  out <- evidence()
  expect_identical(out$lipid_species[4], "MG 18:1")
  expect_identical(out$lipid_source[4], "swisslipids_fallback")
})

test_that("a group row whose key disagrees with its chain sum, and alternative names, get no profile", {
  # Row 5 is an SM group summing to 38:3 but carries a DG key. Since v0.2
  # the key of a group row is looked up, so sl_tier is set, but the
  # Species it gives (DG 36:3) does not match the name, so no profile.
  out <- evidence()
  expect_true(is.na(out$lipid_species[5]))
  expect_false(is.na(out$sl_tier[5]))
  expect_true(is.na(out$lipid_species[8]))
})

test_that("a row found nowhere gets no profile", {
  out <- evidence()
  expect_true(is.na(out$lipid_species[6]))
  expect_false(out$lipid_conflict[6])
})

test_that("disagreeing sources flag a conflict and give no profile", {
  out <- evidence()
  expect_true(out$lipid_conflict[7])
  expect_true(is.na(out$lipid_species[7]))
  expect_true(is.na(out$lipid_level[7]))
})

test_that("row_evidence() must run first", {
  db <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  expect_error(lipid_evidence(make_rows(), db), "row_evidence")
})
