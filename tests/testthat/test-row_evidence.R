test_that("evidence columns are added and inputs are kept", {
  rows <- data.frame(
    row_id        = c("Broad_2022Aug:354", "LEOCC_Metabolon:1"),
    original_name = c("Inosine-15N4 [iSTD]", "bilirubin (E,E)*"),
    inchikey      = c("inchikey", "BPYKTIZUTYGOLE-KDUUSRDASA-N"),
    stringsAsFactors = FALSE
  )
  out <- row_evidence(rows)

  expect_identical(out$row_id, rows$row_id)
  expect_identical(out$inchikey, rows$inchikey)
  expect_identical(out$ik_valid, c(FALSE, TRUE))
  expect_identical(out$ik_block1, c(NA, "BPYKTIZUTYGOLE"))
  expect_identical(out$name_internal_standard, c(TRUE, FALSE))
  expect_identical(out$name_tentative, c(FALSE, TRUE))
  expect_identical(out$name_states_stereo, c(FALSE, TRUE))
})

test_that("a table without an inchikey column still works", {
  out <- row_evidence(data.frame(original_name = "SM 36:1;O2"))
  expect_false(out$ik_valid)
  expect_true(is.na(out$ik_block1))
})

test_that("bad input stops with a clear message", {
  expect_error(row_evidence(list(original_name = "x")), "data frame")
  expect_error(row_evidence(data.frame(name = "x")), "original_name")
  expect_error(
    row_evidence(data.frame(original_name = "x", ik_valid = TRUE)),
    "already has"
  )
})
