test_that("a group takes its most severe pair", {
  # The SM 36:1 group from the review: one broader/narrower pair and two
  # distinct pairs with the SM 38:3 group row.
  pairs <- data.frame(
    group_id = c("SM 36:1;O2", "SM 36:1;O2", "SM 36:1;O2", "Inosine"),
    row_a = c("broad", "broad", "leocc_sm", "a"),
    row_b = c("leocc_sm", "leocc_group", "leocc_group", "b"),
    relationship = c("broader_narrower", "distinct", "distinct", "exact"),
    rule_id = c("R05", "R04", "R04", "R03"),
    stringsAsFactors = FALSE)
  g <- summarise_groups(pairs)
  sm <- g[g$group_id == "SM 36:1;O2", ]
  expect_identical(sm$n_rows, 3L)
  expect_identical(sm$n_pairs, 3L)
  expect_identical(sm$worst_relationship, "distinct")
  expect_identical(sm$worst_rule, "R04")
  expect_identical(sm$n_distinct, 2L)
  expect_identical(sm$n_broader_narrower, 1L)
  expect_identical(g$worst_relationship[g$group_id == "Inosine"], "exact")
})

test_that("unresolved ranks above broader/narrower and below distinct", {
  pairs <- data.frame(group_id = "g", row_a = c("a", "a"), row_b = c("b", "c"),
                      relationship = c("broader_narrower", "unresolved"),
                      rule_id = c("R05", "R11"), stringsAsFactors = FALSE)
  expect_identical(summarise_groups(pairs)$worst_relationship, "unresolved")
})

test_that("an unknown relationship is refused", {
  pairs <- data.frame(group_id = "g", row_a = "a", row_b = "b",
                      relationship = "group_overlap", rule_id = "R08")
  expect_error(summarise_groups(pairs), "Unknown relationship")
})

test_that("same_name ranks below unresolved and above broader/narrower", {
  pairs <- data.frame(group_id = c("g1", "g1", "g2", "g2"), row_a = c("a", "a", "c", "c"),
                      row_b = c("b", "x", "d", "y"),
                      relationship = c("same_name", "unresolved", "same_name", "broader_narrower"),
                      rule_id = c("R01s", "R11", "R01s", "R09"), stringsAsFactors = FALSE)
  g <- summarise_groups(pairs)
  expect_identical(g$worst_relationship, c("unresolved", "same_name"))
  expect_identical(g$n_same_name, c(1L, 1L))
})
