test_that("standard keys are accepted", {
  # Bilirubin, from the LEOCC Metabolon file.
  expect_true(is_valid_inchikey("BPYKTIZUTYGOLE-KDUUSRDASA-N"))
})

test_that("the placeholder text found in the Broad file is rejected", {
  # Broad_2022Aug row 354, "Inosine-15N4 [iSTD]", has the literal text
  # "inchikey" in its InChIKey column.
  expect_false(is_valid_inchikey("inchikey"))
})

test_that("malformed values are rejected, not fixed", {
  expect_false(is_valid_inchikey("InChIKey=BPYKTIZUTYGOLE-KDUUSRDASA-N"))
  expect_false(is_valid_inchikey("bpyktizutygole-kduusrdasa-n"))
  expect_false(is_valid_inchikey("BPYKTIZUTYGOLE-KDUUSRDASA"))
  expect_false(is_valid_inchikey(" BPYKTIZUTYGOLE-KDUUSRDASA-N"))
  expect_false(is_valid_inchikey(""))
})

test_that("missing values give FALSE, never NA", {
  out <- is_valid_inchikey(c(NA, "BPYKTIZUTYGOLE-KDUUSRDASA-N", NA_character_))
  expect_identical(out, c(FALSE, TRUE, FALSE))
})

test_that("output length matches input length", {
  expect_length(is_valid_inchikey(character(0)), 0)
  expect_length(is_valid_inchikey(rep("inchikey", 5)), 5)
})
