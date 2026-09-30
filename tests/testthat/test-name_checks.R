# Every name here is copied from the MetLinkR example files, so the
# tests describe how the checks behave on real data.

test_that("group names are recognised", {
  groups <- c(
    "(2 or 3)-decenoate (10:1n7 or n8)",
    "13-HODE + 9-HODE",
    "2-Hydroxycaproic and or 2-Hydroxyisocaproic acid",
    "10 and or 14-HDoHE",
    "butyrate/isobutyrate (4:0)",
    "sphingomyelin (d18:1/20:2, d18:2/20:1, d16:1/22:2)*",
    "diacylglycerol (16:1/18:2 [2], 16:0/18:3 [1])*",
    "linolenate [alpha or gamma; (18:3n3 or 6)]",
    "aconitate [cis or trans]",
    "(N(1) + N(8))-acetylspermidine",
    "Isobar: fructose 1,6-diphosphate, glucose 1,6-diphosphate, myo-inositol 1,4 or 1,3-diphosphate"
  )
  expect_true(all(is_group_name(groups)))
})

test_that("a slash between lipid chains is not a group", {
  single <- c(
    "stearoyl sphingomyelin (d18:1/18:0)",
    "oleoyl-linoleoyl-glycerol (18:1/18:2) [1]",
    "SM 36:1;O2",
    "asparagine"
  )
  expect_false(any(is_group_name(single)))
})

test_that("top-level 'or' between full names is kept apart from groups", {
  alt <- c(
    "Sarcosine or Sarcosine (N-Methylglycine)",
    "1-palmitoyl-GPC (16:0) or 1-palmitoylglycerophosphocholine (16:0)",
    "branched-chain, straight-chain, or cyclopropyl 12:1 fatty acid*"
  )
  expect_true(all(has_alternative_names(alt)))
  expect_false(any(is_group_name(alt)))
})

test_that("'and or', bracketed 'or' and locant 'or' are not alternatives", {
  not_alt <- c(
    "Citric and or Isocitric acid",
    "(12 or 13)-methylmyristate (a15:0 or i15:0)",
    "Isobar: fructose 1,6-diphosphate, glucose 1,6-diphosphate, myo-inositol 1,4 or 1,3-diphosphate",
    "asparagine"
  )
  expect_false(any(has_alternative_names(not_alt)))
})

test_that("tentative marks are found", {
  expect_true(is_tentative_name("andro steroid monosulfate C19H28O6S (1)*"))
  expect_true(is_tentative_name("hydroxypalmitoyl sphingomyelin (d18:1/16:0(OH))**"))
  expect_false(is_tentative_name("androsterone sulfate"))
})

test_that("isomer labels are found only at the end, after a space", {
  expect_true(all(has_isomer_label(c(
    "2-methoxyhydroquinone sulfate (1)",
    "oleoyl-linoleoyl-glycerol (18:1/18:2) [2]",
    "andro steroid monosulfate C19H28O6S (1)*"
  ))))
  # Residue positions in peptide names are not isomer labels.
  expect_false(any(has_isomer_label(c(
    "bradykinin, des-arg(9)",
    "Fibrinopeptide A, des-ala(1)*",
    "stearoyl sphingomyelin (d18:1/18:0)"
  ))))
})

test_that("internal standards are found", {
  expect_true(all(is_internal_standard(c(
    "Inosine-15N4 [iSTD]",
    "Glycocholate-d4 [iSTD]",
    "Thymine-d4 [iSTD]",
    "15S-15-methyl PGE2 [iSTD]"
  ))))
})

test_that("formulas and sphingoid bases are not mistaken for isotope labels", {
  expect_false(any(is_internal_standard(c(
    "glucuronide of C10H18O2 (1)*",
    "stearoyl sphingomyelin (d18:1/18:0)",
    "N-palmitoyl-sphingadienine (d18:2/16:0)*",
    "inosine"
  ))))
})

test_that("stated stereo is found", {
  expect_true(all(states_stereo(c(
    "bilirubin (E,E)*",
    "bilirubin (Z,Z)",
    "cis-4-decenoylcarnitine (C10:1)",
    "trans-4-hydroxyproline",
    "(S)-3-hydroxybutyrylcarnitine",
    "(R)-3-hydroxybutyrylcarnitine",
    "Leucine (L-leucine)",
    "O-sulfo-L-tyrosine"
  ))))
})

test_that("positional alpha and beta, racemic marks and plain names do not count", {
  expect_false(any(states_stereo(c(
    "alpha-hydroxyisovalerate",
    "alpha-ketoglutarate",
    "2-hydroxystearate",
    "2-Hydroxyglutaric acid",
    "DL-lactate",
    "HDL-cholesterol"
  ))))
})

test_that("missing names give FALSE everywhere", {
  expect_identical(is_group_name(NA), FALSE)
  expect_identical(has_alternative_names(NA), FALSE)
  expect_identical(is_tentative_name(NA), FALSE)
  expect_identical(has_isomer_label(NA), FALSE)
  expect_identical(is_internal_standard(NA), FALSE)
  expect_identical(states_stereo(NA), FALSE)
})
