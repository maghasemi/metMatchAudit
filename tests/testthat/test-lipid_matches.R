# Evidence is built directly, so these tests do not depend on SwissLipids
# or rgoslin. Names and Species follow reviewed proposals from the
# example data (notes/v02_proposal_review.csv).
match_evidence <- function() {
  rows <- data.frame(
    row_id = c("L1", "L2", "L3", "L4", "B1", "B2", "B3"),
    file   = c("LEOCC", "LEOCC", "LEOCC", "LEOCC", "Broad", "Broad", "Broad"),
    original_name = c("1-myristoyl-2-arachidonoyl-GPC (14:0/20:4)*",
                      "1-palmitoleoyl-2-linolenoyl-GPC (16:1/18:3)*",
                      "1-(1-enyl-palmitoyl)-2-oleoyl-GPC (P-16:0/18:1)*",
                      "stearoyl sphingomyelin (d18:1/18:0)",
                      "PC 34:4", "PC O-34:2", "SM 36:1;O2"),
    inchikey = NA_character_,
    lipid_species = c("PC 34:4", "PC 34:4", "PC O-34:2", "SM 36:1;O2",
                      "PC 34:4", "PC O-34:2", "SM 36:1;O2"),
    lipid_molecular = NA_character_,
    lipid_level = c("structure", "structure", "structure", "structure",
                    "species", "species", "species"),
    group_id = c(NA, NA, NA, "SM 18:1;O2/18:0", NA, NA, "SM 18:1;O2/18:0"),
    stringsAsFactors = FALSE)
  row_evidence(rows)
}

test_that("matches are proposed across files with the right relationship", {
  m <- propose_lipid_matches(match_evidence())
  expect_equal(nrow(m), 4)
  expect_true(all(m$relationship == "broader_narrower"))
  expect_true(all(m$file_a != m$file_b))
})

test_that("one Broad Species against two Metabolon structures is one_to_many", {
  m <- propose_lipid_matches(match_evidence())
  pc <- m[m$lipid_species == "PC 34:4", ]
  expect_equal(nrow(pc), 2)
  expect_true(all(pc$shape == "one_to_many"))
  expect_true(all(pc$partners_b == 2))
})

test_that("a plasmalogen matches the O notation Species one to one", {
  m <- propose_lipid_matches(match_evidence())
  expect_identical(m$shape[m$lipid_species == "PC O-34:2"], "one_to_one")
})

test_that("existing groups are reported, new proposals are not marked grouped", {
  m <- propose_lipid_matches(match_evidence())
  expect_true(m$already_grouped[m$lipid_species == "SM 36:1;O2"])
  expect_false(any(m$already_grouped[m$lipid_species != "SM 36:1;O2"]))
})

test_that("no group column is needed", {
  ev <- match_evidence()
  ev$group_id <- NULL
  m <- propose_lipid_matches(ev, group_col = NULL)
  expect_equal(nrow(m), 4)
  expect_false("already_grouped" %in% names(m))
})

test_that("rows of the same file are never paired", {
  ev <- match_evidence()
  ev$file <- "one_file"
  expect_equal(nrow(propose_lipid_matches(ev)), 0)
})
