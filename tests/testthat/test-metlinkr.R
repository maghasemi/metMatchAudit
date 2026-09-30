mini_files <- function() {
  data.frame(
    path = test_path("fixtures", "metlinkr",
                     c("LEOCC_mini_metLinkR.csv", "Broad_mini_metLinkR.csv")),
    name_col = c("CHEMICAL_NAME", "Metabolite"),
    inchikey_col = c("INCHIKEY", "inChIKey"),
    stringsAsFactors = FALSE)
}

test_that("MetLinkR files are read into one row table", {
  rows <- read_metlinkr(mini_files())
  expect_equal(nrow(rows), 9)
  expect_identical(rows$row_id[1], "LEOCC_mini:1")
  expect_identical(rows$row_id[6], "Broad_mini:1")
  expect_identical(rows$group_id[1], "SM 18:1;O2/18:0")
})

test_that("NA, empty and '-' standardized names mean no group", {
  rows <- read_metlinkr(mini_files())
  expect_true(is.na(rows$group_id[rows$original_name == "cis-4-decenoylcarnitine (C10:1)"]))
  expect_true(is.na(rows$group_id[rows$original_name == "Glycine"]))
})

test_that("a file without an InChIKey column is accepted", {
  f <- mini_files()
  f$inchikey_col[2] <- NA
  rows <- read_metlinkr(f)
  expect_true(all(is.na(rows$inchikey[rows$file == "Broad_mini"])))
})

test_that("a wrong column name stops with the file and column named", {
  f <- mini_files()
  f$name_col[1] <- "BIOCHEMICAL"
  expect_error(read_metlinkr(f), "BIOCHEMICAL")
})

test_that("the full audit reproduces the reviewed results on these rows", {
  db <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  a <- audit_metlinkr(mini_files(), db)
  expect_s3_class(a, "metmatch_audit")

  rule_for <- function(x, y) {
    ids <- a$rows$row_id[match(c(x, y), a$rows$original_name)]
    p <- a$pairs[(a$pairs$row_a == ids[1] & a$pairs$row_b == ids[2]) |
                 (a$pairs$row_a == ids[2] & a$pairs$row_b == ids[1]), ]
    p$rule_id
  }
  expect_identical(rule_for("stearoyl sphingomyelin (d18:1/18:0)", "SM 36:1;O2"), "R05")
  expect_identical(rule_for("inosine", "Inosine"), "R03")
  expect_identical(rule_for("Inosine", "Inosine-15N4 [iSTD]"), "R02")
  expect_identical(rule_for("oleoyl-linoleoyl-glycerol (18:1/18:2) [1]",
                            "oleoyl-linoleoyl-glycerol (18:1/18:2) [2]"), "R12")

  expect_identical(a$row_flags$check, "C01")
  expect_identical(nrow(a$groups), 3L)
  expect_identical(a$info$swisslipids_md5, db$md5)
})

test_that("write_audit writes the five files", {
  db <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  a <- audit_metlinkr(mini_files(), db)
  dir <- tempfile()
  paths <- write_audit(a, dir)
  expect_true(all(file.exists(paths)))
  expect_equal(nrow(utils::read.csv(file.path(dir, "pairs.csv"))), nrow(a$pairs))
})

test_that("print shows a summary without error", {
  db <- load_swisslipids(test_path("fixtures", "swisslipids_mini.tsv"))
  a <- audit_metlinkr(mini_files(), db)
  expect_output(print(a), "metmatch_audit")
})
