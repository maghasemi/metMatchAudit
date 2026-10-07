# LMSD lookup on a small SDF written here, in the layout of the LMSD
# download: a molfile block, then "> <FIELD>" lines, then "$$$$".

lmsd_fixture <- function() {
  rec <- function(id, key, abbr) {
    c(id, "  LIPDMAPS", "", "  0  0  0  0  0  0  0  0  0  0999 V2000", "M  END",
      "> <LM_ID>", id, "", "> <INCHI_KEY>", key, "",
      if (!is.na(abbr)) c("> <ABBREVIATION>", abbr, ""), "$$$$")
  }
  f <- tempfile(fileext = ".sdf")
  writeLines(c(
    rec("LMFA07070001", "NNCBVXBBLABOCB-SEYXRHQNSA-N", "CAR 14:1"),
    # Same connectivity, other stereo, same Species.
    rec("LMFA07070002", "NNCBVXBBLABOCB-UHFFFAOYSA-N", "CAR 14:1"),
    # Two records sharing block1 with different Species.
    rec("LMGP01010001", "AAAAAAAAAAAAAA-UHFFFAOYSA-N", "PC 34:1"),
    rec("LMGP01010002", "AAAAAAAAAAAAAA-BBBBBBBBBB-N", "PC 34:2"),
    # No abbreviation: not used.
    rec("LMST01010001", "CCCCCCCCCCCCCC-UHFFFAOYSA-N", NA)
  ), f)
  load_lmsd(f)
}

test_that("load_lmsd keeps records with a key and an abbreviation", {
  db <- lmsd_fixture()
  expect_s3_class(db, "lmsd_db")
  expect_identical(db$n_records, 5L)
  expect_identical(db$n_used, 4L)
  expect_error(load_lmsd(tempfile()), "Cannot find")
})

test_that("lmsd_lookup uses the first tier that finds anything", {
  out <- lmsd_lookup(c("NNCBVXBBLABOCB-SEYXRHQNSA-N",   # full key
                       "NNCBVXBBLABOCB-SEYXRHQNSA-M",   # first 25 characters
                       "NNCBVXBBLABOCB-ZZZZZZZZZZ-N",   # first 14
                       "DDDDDDDDDDDDDD-UHFFFAOYSA-N",   # nothing
                       NA), lmsd_fixture())
  expect_identical(out$lm_tier, c(1L, 2L, 3L, NA, NA))
  expect_identical(out$lm_species[1:3], rep("CAR 14:1", 3))
  expect_true(is.na(out$lm_species[4]))
})

test_that("hits that disagree on Species give no Species", {
  out <- lmsd_lookup("AAAAAAAAAAAAAA-ZZZZZZZZZZ-N", lmsd_fixture())
  expect_identical(out$lm_n_hits, 2L)
  expect_true(is.na(out$lm_species))
})

test_that("lipid_evidence uses LMSD only where SwissLipids and the name give nothing", {
  sl <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  rows <- data.frame(
    row_id = c("a", "b"),
    original_name = c("myristoleoylcarnitine (C14:1)*", "SM 36:1;O2"),
    inchikey = c("NNCBVXBBLABOCB-SEYXRHQNSA-N", "LKQLRGMMMAHREN-YJFXYUILSA-N"),
    stringsAsFactors = FALSE)
  ev <- row_evidence(rows)
  without <- lipid_evidence(ev, sl)
  with    <- lipid_evidence(ev, sl, lmsd = lmsd_fixture())
  expect_true(is.na(without$lipid_species[1]))
  expect_identical(with$lipid_species[1], "CAR 14:1")
  expect_identical(with$lipid_source[1], "lmsd")
  expect_identical(with$lipid_level[1], "structure")
  # The SM row keeps the profile it had.
  expect_identical(with[2, c("lipid_species", "lipid_source", "lipid_level")],
                   without[2, c("lipid_species", "lipid_source", "lipid_level")])
})
