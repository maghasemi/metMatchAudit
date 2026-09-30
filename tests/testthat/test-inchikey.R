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

test_that("keys split into the right parts", {
  p <- inchikey_parts("BPYKTIZUTYGOLE-KDUUSRDASA-N")
  expect_identical(p$block1, "BPYKTIZUTYGOLE")
  expect_identical(p$block2, "KDUUSRDASA")
  expect_identical(p$protonation, "N")
  expect_false(p$stereo_unspecified)
})

test_that("UHFFFAOYSA marks a key without stereo", {
  # 9-hydroxystearate in the LEOCC file.
  p <- inchikey_parts("RKHXDCVAPIMDMG-UHFFFAOYSA-N")
  expect_true(p$stereo_unspecified)
})

test_that("invalid keys give NA parts, not partial strings", {
  p <- inchikey_parts(c("inchikey", NA, "BPYKTIZUTYGOLE-KDUUSRDASA-N"))
  expect_identical(p$block1, c(NA, NA, "BPYKTIZUTYGOLE"))
  expect_identical(p$stereo_unspecified, c(NA, NA, FALSE))
  expect_equal(nrow(p), 3)
})

test_that("keys that differ only in stereo share block1", {
  # 2-hydroxystearic acid from the review: 2S in one file, 2R in the
  # other (LEOCC and Broad).
  p <- inchikey_parts(c("KIHBGTRZFAVZRV-QGZVFWFLSA-N",
                        "KIHBGTRZFAVZRV-KRWDZBQOSA-N"))
  expect_identical(p$block1[1], p$block1[2])
  expect_false(identical(p$block2[1], p$block2[2]))
})
