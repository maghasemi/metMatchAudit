flag_input <- function() {
  rows <- data.frame(
    row_id = c("Broad_2022Aug:354", "LEOCC:inosine", "Broad:nokey", "made-up:conflict"),
    original_name = c("Inosine-15N4 [iSTD]", "inosine", "Sarcosine", "DG 18:1_18:2"),
    inchikey = c("inchikey", "UGQMRVRMYYASKQ-KQYNXXCUSA-N", NA, "LKQLRGMMMAHREN-YJFXYUILSA-N"),
    stringsAsFactors = FALSE)
  ev <- row_evidence(rows)
  ev$lipid_conflict <- c(FALSE, FALSE, FALSE, TRUE)
  ev$sl_species_abbrev <- c(NA, NA, NA, "SM(d36:1)")
  ev
}

test_that("C01 flags the placeholder key from the Broad file", {
  f <- row_flags(flag_input())
  c01 <- f[f$check == "C01", ]
  expect_identical(c01$row_id, "Broad_2022Aug:354")
  expect_match(c01$detail, "inchikey")
})

test_that("a missing key is not flagged", {
  f <- row_flags(flag_input())
  expect_false("Broad:nokey" %in% f$row_id)
})

test_that("C03 flags a lipid source conflict", {
  f <- row_flags(flag_input())
  c03 <- f[f$check == "C03", ]
  expect_identical(c03$row_id, "made-up:conflict")
  expect_match(c03$detail, "SM\\(d36:1\\)")
})

test_that("clean rows give an empty table with the right columns", {
  ev <- flag_input()[2, ]
  f <- row_flags(ev)
  expect_equal(nrow(f), 0)
  expect_named(f, c("row_id", "check", "flag", "detail"))
})

test_that("C04 separates connectivity from stereo differences", {
  rows <- data.frame(
    row_id = c("Broad_2022Aug:384", "Broad_2022Aug:504", "ok"),
    original_name = c("Xylose", "N-Formylmethionine", "inosine"),
    inchikey = c("SRBFZHDQGSBBOR-IOVATXLUSA-N", "PYUSHNKNPOHWEZ-RXMQYKEDSA-N",
                 "UGQMRVRMYYASKQ-KQYNXXCUSA-N"),
    stringsAsFactors = FALSE)
  ev <- row_evidence(rows)
  ev$lipid_conflict <- FALSE
  ev$sl_species_abbrev <- NA
  ev$st_key <- c("PYMYPHUHKUWMLA-VPENINKCSA-N", "PYUSHNKNPOHWEZ-YFKPBYRVSA-N",
                 "UGQMRVRMYYASKQ-KQYNXXCUSA-N")
  ev$st_key_agrees <- c(FALSE, FALSE, TRUE)
  ev$st_block1_differs <- c(TRUE, FALSE, FALSE)
  f <- row_flags(ev)
  expect_identical(f$flag[f$row_id == "Broad_2022Aug:384"], "structure_key_connectivity")
  expect_identical(f$flag[f$row_id == "Broad_2022Aug:504"], "structure_key_stereo")
  expect_false("ok" %in% f$row_id)
})
