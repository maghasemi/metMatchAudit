test_that("/m1 flips the signs of /t", {
  s <- inchi_stereo(inchi_cholesterol_broad)
  # Written as 21+ with /m1, so the real sign is minus.
  expect_identical(unname(s$centres["t21"]), "-")
  expect_identical(s$s, "1")
})

test_that("undefined centres are dropped", {
  s <- inchi_stereo(inchi_cholesterol_leocc)
  expect_identical(names(s$centres), "t21")
})

test_that("cholesterol with one defined centre is a subset of full cholesterol", {
  expect_identical(compare_stereo(inchi_cholesterol_leocc, inchi_cholesterol_broad),
                   "a_subset_of_b")
  expect_identical(compare_stereo(inchi_cholesterol_broad, inchi_cholesterol_leocc),
                   "b_subset_of_a")
})

test_that("the LEOCC taurocholate structure conflicts with the Broad one (C9)", {
  expect_identical(compare_stereo(inchi_taurocholate_leocc, inchi_taurocholate_broad),
                   "conflict")
})

test_that("double bond geometry is compared too (adrenic acid, all-trans vs all-cis)", {
  expect_identical(compare_stereo(inchi_adrenate_leocc, inchi_adrenate_broad), "conflict")
  expect_length(inchi_stereo(inchi_adrenate_broad)$centres, 4)
})

test_that("identical stereo is same, and relative stereo is not compared", {
  expect_identical(compare_stereo(inchi_cholesterol_broad, inchi_cholesterol_broad), "same")
  rel <- sub("/s1$", "/s2", inchi_cholesterol_broad)
  expect_identical(compare_stereo(rel, inchi_cholesterol_broad), "relative")
})

test_that("InChIs without stereo, and non-InChI input", {
  no_st <- "InChI=1S/C4H10O4/c5-1-3(7)4(8)2-6/h3-8H,1-2H2"
  expect_identical(compare_stereo(no_st, no_st), "no_stereo")
  expect_identical(compare_stereo(no_st, inchi_xylose_broad), "a_subset_of_b")
  expect_length(inchi_stereo(NA)$centres, 0)
  expect_length(inchi_stereo("not an inchi")$centres, 0)
})

test_that("layers after a fixed-H layer are ignored", {
  x <- paste0(inchi_cholesterol_leocc, "/f/h/t19-")
  expect_identical(names(inchi_stereo(x)$centres), "t21")
})
