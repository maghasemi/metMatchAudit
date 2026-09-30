fixture <- function() test_path("fixtures", "swisslipids_mini.tsv")

test_that("the table loads and records its md5", {
  db <- load_swisslipids(fixture())
  expect_s3_class(db, "swisslipids_db")
  expect_equal(db$n_rows, 25)
  expect_identical(db$md5, unname(tools::md5sum(fixture())))
})

test_that("a wrong md5 stops loading", {
  expect_error(load_swisslipids(fixture(), expected_md5 = "0000"), "md5")
})

test_that("an unexpected header stops loading", {
  tmp <- tempfile(fileext = ".tsv")
  lines <- readLines(fixture())
  lines[1] <- sub("Parent", "Parents", lines[1], fixed = TRUE)
  writeLines(lines, tmp)
  expect_error(load_swisslipids(tmp), "header")
})

test_that("a full key match is tier 1 and climbs to Species and Molecular", {
  db <- load_swisslipids(fixture())
  out <- swisslipids_lookup("BLZVZPYMHLXLHG-RQOIEFAZSA-N", db)
  expect_identical(out$sl_tier, 1L)
  expect_identical(out$sl_species, "DG(36:3)")
  expect_identical(out$sl_molecular, "DG(18:1_18:2)")
  expect_false(out$sl_fallback)
})

test_that("a different protonation letter matches at tier 2", {
  db <- load_swisslipids(fixture())
  out <- swisslipids_lookup("GREDRAMJRDQWEJ-LTEAFHAISA-N", db)
  expect_identical(out$sl_tier, 2L)
  expect_identical(out$sl_species, "DG(36:3)")
})

test_that("a different stereo block matches at tier 3", {
  db <- load_swisslipids(fixture())
  out <- swisslipids_lookup("BLZVZPYMHLXLHG-ZZZZZZZZZZ-N", db)
  expect_identical(out$sl_tier, 3L)
  expect_identical(out$sl_species, "DG(36:3)")
})

test_that("an entry without a hierarchy uses the connectivity fallback", {
  # Seen with 1-oleoylglycerol in the feasibility run.
  db <- load_swisslipids(fixture())
  out <- swisslipids_lookup("RZRNAYUHWVFMIP-KTKRTIGZSA-N", db)
  expect_identical(out$sl_tier, 1L)
  expect_identical(out$sl_species, "MG(18:1)")
  expect_true(out$sl_fallback)
})

test_that("the fallback refuses when siblings name two Species", {
  db <- load_swisslipids(fixture())
  out <- swisslipids_lookup("AAAAAAAAAAAAAA-KTKRTIGZSA-N", db)
  expect_identical(out$sl_tier, 1L)
  expect_true(is.na(out$sl_species))
  expect_identical(out$sl_n_species, 0L)
  expect_false(out$sl_fallback)
})

test_that("a loop in the parent links does not hang", {
  db <- load_swisslipids(fixture())
  out <- swisslipids_lookup("LOOPLOOPLOOPLO-KTKRTIGZSA-N", db)
  expect_identical(out$sl_tier, 1L)
  expect_true(is.na(out$sl_species))
})

test_that("unknown and invalid keys give an empty row", {
  db <- load_swisslipids(fixture())
  out <- swisslipids_lookup(c("OQWOHRPOYAVIOK-KTKRTIGZSA-N", "inchikey", NA), db)
  expect_equal(nrow(out), 3)
  expect_true(all(is.na(out$sl_tier)))
  expect_identical(out$sl_n_hits, c(0L, 0L, 0L))
})

test_that("lookup refuses anything that is not a loaded table", {
  expect_error(swisslipids_lookup("X", list()), "load_swisslipids")
})
