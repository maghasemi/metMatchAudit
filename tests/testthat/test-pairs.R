test_that("a group of three gives three pairs, in input order", {
  rows <- data.frame(row_id = c("r1", "r2", "r3", "r4"),
                     group_id = c("SM 36:1", "SM 36:1", "SM 36:1", "Xylose"))
  p <- make_pairs(rows)
  expect_equal(nrow(p), 3)
  expect_identical(p$row_a, c("r1", "r1", "r2"))
  expect_identical(p$row_b, c("r2", "r3", "r3"))
})

test_that("single rows and missing groups give no pairs", {
  rows <- data.frame(row_id = c("r1", "r2", "r3"), group_id = c("A", NA, ""))
  expect_equal(nrow(make_pairs(rows)), 0)
})

test_that("duplicate row IDs are refused", {
  expect_error(make_pairs(data.frame(row_id = c("x", "x"), group_id = "A")),
               "unique")
})
