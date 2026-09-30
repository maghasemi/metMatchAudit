test_that("members of real Metabolon group names share one chain sum", {
  expect_identical(
    group_chain_sum(c("sphingomyelin (d18:1/20:2, d18:2/20:1, d16:1/22:2)*",
                      "sphingomyelin (d18:1/25:0, d19:0/24:1, d20:1/23:0, d19:1/24:0)*",
                      "diacylglycerol (16:1/18:2 [2], 16:0/18:3 [1])*")),
    c("38:3", "43:1", "34:3"))
})

test_that("members that disagree, unreadable members and plain names give NA", {
  expect_true(is.na(group_chain_sum("sphingomyelin (d18:1/20:2, d18:1/20:1)*")))
  expect_true(is.na(group_chain_sum("sphingomyelin (d18:1/16:0(OH), d18:1/16:0)*")))
  expect_true(is.na(group_chain_sum("stearoyl sphingomyelin (d18:1/18:0)")))
  expect_true(is.na(group_chain_sum(NA)))
})

group_rows <- function() {
  data.frame(
    row_id = c("g38", "g34", "bad"),
    original_name = c("sphingomyelin (d18:1/20:2, d18:2/20:1, d16:1/22:2)*",
                      "sphingomyelin (d18:2/16:0, d18:1/16:1)*",
                      # Name says 38:3, but the key is the one of an SM 36:1.
                      "sphingomyelin (d18:1/20:2, d18:2/20:1)*"),
    inchikey = c("YVTVVJRJLRHYSJ-COUPUTHFSA-N", "IWHLTJGDYLJQSH-FLPFCPAFSA-N",
                 "LKQLRGMMMAHREN-YJFXYUILSA-N"),
    stringsAsFactors = FALSE)
}

test_that("a group row gets a group-level profile when both routes agree", {
  db <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  ev <- lipid_evidence(row_evidence(group_rows()), db)
  expect_identical(ev$lipid_species[1:2], c("SM 38:3;O2", "SM 34:2;O2"))
  expect_identical(ev$lipid_level[1:2], c("group", "group"))
  expect_identical(ev$lipid_source[1:2], c("group_chain_sum", "group_chain_sum"))
})

test_that("a group row whose key points to another Species gets no profile", {
  db <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  ev <- lipid_evidence(row_evidence(group_rows()), db)
  expect_true(is.na(ev$lipid_species[3]))
})
