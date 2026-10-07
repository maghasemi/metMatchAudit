# One test per rule, built from pairs seen in the manual review of the
# MetLinkR example output (notes/review.csv) wherever possible. Names
# and InChIKeys are copied from the source files.

test_that("R02: an internal standard is never the same as the analyte", {
  # Inosine-15N4 [iSTD] was grouped with inosine. Its InChIKey field
  # holds the text "inchikey", so R02 must run before R01.
  r <- pair_result(list(name = "inosine", inchikey = "UGQMRVRMYYASKQ-KQYNXXCUSA-N"),
                   list(name = "Inosine-15N4 [iSTD]", inchikey = "inchikey"))
  expect_identical(r$rule_id, "R02")
  expect_identical(r$relationship, "distinct")
})

test_that("R01: no key and no lipid profile leaves the pair unresolved", {
  r <- pair_result(list(name = "Sarcosine"),
                   list(name = "sarcosine", inchikey = "FSYKKLYZXJSNPZ-UHFFFAOYSA-N"))
  expect_identical(r$rule_id, "R01")
  expect_identical(r$relationship, "unresolved")
})

test_that("R05 runs before R03: same key, but a Species against a structure", {
  # Broad reports "SM 36:1;O2" and attached the same InChIKey as the
  # Metabolon row "stearoyl sphingomyelin (d18:1/18:0)". Manual label: D.
  key <- "LKQLRGMMMAHREN-YJFXYUILSA-N"
  r <- pair_result(
    list(name = "SM 36:1;O2", inchikey = key,
         lipid_species = "SM 36:1;O2", lipid_level = "species"),
    list(name = "stearoyl sphingomyelin (d18:1/18:0)", inchikey = key,
         lipid_species = "SM 36:1;O2", lipid_molecular = "SM 18:1;O2/18:0",
         lipid_level = "structure"))
  expect_identical(r$rule_id, "R05")
  expect_identical(r$relationship, "broader_narrower")
})

test_that("R04: different lipid Species are different compounds", {
  r <- pair_result(
    list(name = "SM 36:1;O2", lipid_species = "SM 36:1;O2", lipid_level = "species"),
    list(name = "SM 38:1;O2", lipid_species = "SM 38:1;O2", lipid_level = "species"))
  expect_identical(r$rule_id, "R04")
  expect_identical(r$relationship, "distinct")
})

test_that("R04h: HexCer against GlcCer with the same chains stays unresolved", {
  # The four hierarchy vs hit disagreements in the feasibility run.
  r <- pair_result(
    list(name = "glycosyl-N-palmitoyl-sphingosine (d18:1/16:0)",
         inchikey = "VJLLLMIZEJJZTE-NNTBDIJYSA-N",
         lipid_species = "HexCer 34:1;O2", lipid_level = "structure"),
    list(name = "GlcCer 34:1;O2",
         lipid_species = "GlcCer 34:1;O2", lipid_level = "species"))
  expect_identical(r$rule_id, "R04h")
  expect_identical(r$relationship, "unresolved")
})

test_that("R04h does not hide a real class difference such as PC against PE", {
  r <- pair_result(
    list(name = "PC 34:1", lipid_species = "PC 34:1", lipid_level = "species"),
    list(name = "PE 34:1", lipid_species = "PE 34:1", lipid_level = "species"))
  expect_identical(r$rule_id, "R04")
})

test_that("R06: same Species, different molecular species", {
  r <- pair_result(
    list(name = "PC 16:0_18:1", lipid_species = "PC 34:1",
         lipid_molecular = "PC 16:0_18:1", lipid_level = "molecular"),
    list(name = "PC 18:0_16:1", lipid_species = "PC 34:1",
         lipid_molecular = "PC 16:1_18:0", lipid_level = "molecular"))
  expect_identical(r$rule_id, "R06")
  expect_identical(r$relationship, "distinct")
})

test_that("R06b: the same Species reported by two vendors is exact", {
  r <- pair_result(
    list(name = "SM 36:1;O2", lipid_species = "SM 36:1;O2", lipid_level = "species"),
    list(name = "SM(d36:1)", lipid_species = "SM 36:1;O2", lipid_level = "species"))
  expect_identical(r$rule_id, "R06b")
  expect_identical(r$relationship, "exact")
})

test_that("R03: identical keys on non-lipid rows are exact", {
  key <- "UGQMRVRMYYASKQ-KQYNXXCUSA-N"
  r <- pair_result(list(name = "inosine", inchikey = key),
                   list(name = "Inosine", inchikey = key))
  expect_identical(r$rule_id, "R03")
  expect_identical(r$relationship, "exact")
})

test_that("R08: a group name leaves the pair unresolved", {
  r <- pair_result(
    list(name = "(2 or 3)-decenoate (10:1n7 or n8)", inchikey = "WXBXVVIUZANZAU-CMDGGOBGSA-N"),
    list(name = "3-Decenoic acid", inchikey = "CPVUNKGURQKKKX-BQYQJAHWSA-N"))
  expect_identical(r$rule_id, "R08")
  expect_identical(r$relationship, "unresolved")
})

test_that("R08b: alternative names leave the pair unresolved", {
  key <- "FSYKKLYZXJSNPZ-UHFFFAOYSA-N"
  r <- pair_result(list(name = "Sarcosine or Sarcosine (N-Methylglycine)", inchikey = key),
                   list(name = "sarcosine", inchikey = "FSYKKLYZXJSNPZ-UHFFFAOYSA-M"))
  expect_identical(r$rule_id, "R08b")
})

test_that("R09: a key with stereo against one without is broader/narrower", {
  # 2-hydroxyglutarate: no stereo in the LEOCC key, 2S (L) in Broad.
  r <- pair_result(list(name = "2-hydroxyglutarate", inchikey = "HWXBTNAVRSUOJR-UHFFFAOYSA-N"),
                   list(name = "2-Hydroxyglutaric acid", inchikey = "HWXBTNAVRSUOJR-VKHMYHEASA-N"))
  expect_identical(r$rule_id, "R09")
  expect_identical(r$relationship, "broader_narrower")
})

test_that("R10: stereo differences stated in the names", {
  r <- pair_result(list(name = "bilirubin (Z,Z)", inchikey = "BPYKTIZUTYGOLE-KDUUSRDASA-N"),
                   list(name = "bilirubin (E,E)*", inchikey = "BPYKTIZUTYGOLE-VVCLLGATSA-N"))
  expect_identical(r$rule_id, "R10")
  expect_identical(r$relationship, "stereo_variant")
})

test_that("R11: stereo only in the keys stays unresolved", {
  # 2-hydroxystearate as 2S in one file and 2R in the other; neither
  # name states stereo.
  r <- pair_result(list(name = "2-hydroxystearate", inchikey = "KIHBGTRZFAVZRV-QGZVFWFLSA-N"),
                   list(name = "2-Hydroxystearic acid", inchikey = "KIHBGTRZFAVZRV-KRWDZBQOSA-N"))
  expect_identical(r$rule_id, "R11")
  expect_identical(r$relationship, "unresolved")
})

test_that("R12: different connectivity is distinct", {
  # 2-AG and 1-AG, manual label C.
  r <- pair_result(list(name = "2-arachidonoylglycerol (20:4)", inchikey = "RCRCTBLIHCHWDZ-DOFZRALJSA-N"),
                   list(name = "1-arachidonylglycerol (20:4)", inchikey = "DCPCOKIYJYGMDN-DOFZRALJSA-N"))
  expect_identical(r$rule_id, "R12")
  expect_identical(r$relationship, "distinct")
})

test_that("R12 on two structure-level lipids with the same molecular species", {
  # oleoyl-linoleoyl-glycerol [1] and [2]: both DG 36:3 and DG 18:1_18:2,
  # but a 1,2- against a 1,3-diacylglycerol. Manual label C.
  r <- pair_result(
    list(name = "oleoyl-linoleoyl-glycerol (18:1/18:2) [1]", inchikey = "BLZVZPYMHLXLHG-RQOIEFAZSA-N",
         lipid_species = "DG 36:3", lipid_molecular = "DG 18:1_18:2", lipid_level = "structure"),
    list(name = "oleoyl-linoleoyl-glycerol (18:1/18:2) [2]", inchikey = "GREDRAMJRDQWEJ-LTEAFHAISA-N",
         lipid_species = "DG 36:3", lipid_molecular = "DG 18:1_18:2", lipid_level = "structure"))
  expect_identical(r$rule_id, "R12")
})

test_that("R12n: xylose ring and open-chain forms are no longer called distinct", {
  # Until v0.4.1 this pair came out distinct through R12 (a known limit).
  # The names are the same, so R12n now leaves it unresolved.
  r <- pair_result(list(name = "xylose", inchikey = "PYMYPHUHKUWMLA-VPENINKCSA-N"),
                   list(name = "Xylose", inchikey = "SRBFZHDQGSBBOR-IOVATXLUSA-N"))
  expect_identical(r$rule_id, "R12n")
  expect_identical(r$relationship, "unresolved")
})

test_that("R12n: same name with conflicting keys from the COMETS test", {
  # fructose, one key for the furanose ring and one for the open chain.
  r <- pair_result(list(name = "fructose", inchikey = "RFSUNEUAIZKAJO-ARQDHWQXSA-N"),
                   list(name = "fructose", inchikey = "BJHIKXHVCXFQLS-UYFOZJQFSA-N"))
  expect_identical(r$rule_id, "R12n")
  expect_identical(r$relationship, "unresolved")
})

test_that("R12n ignores a trailing star but not an isomer label", {
  r <- pair_result(list(name = "4-vinylphenol sulfate", inchikey = "FUGYGGDSWSUORM-UHFFFAOYSA-N"),
                   list(name = "4-Vinylphenol sulfate*", inchikey = "IETVQHUKTKKBFF-UHFFFAOYSA-N"))
  expect_identical(r$rule_id, "R12n")
  r <- pair_result(list(name = "2-methoxyhydroquinone sulfate (1)", inchikey = "FUGYGGDSWSUORM-UHFFFAOYSA-N"),
                   list(name = "2-methoxyhydroquinone sulfate (2)", inchikey = "IETVQHUKTKKBFF-UHFFFAOYSA-N"))
  expect_identical(r$rule_id, "R12")
})

test_that("R03c: one stereo-free key for scyllo- and myo-inositol is not exact", {
  key <- "CDAISMWEOUEBRE-UHFFFAOYSA-N"
  r <- pair_result(list(name = "scyllo-inositol", inchikey = key),
                   list(name = "myo-inositol", inchikey = key))
  expect_identical(r$rule_id, "R03c")
  expect_identical(r$relationship, "unresolved")
})

test_that("R03c: one key for n3 and n6 DPA is not exact", {
  key <- "AVKOENOBFIYBSA-WMPRHZDHSA-N"
  r <- pair_result(list(name = "docosapentaenoate (n3 DPA; 22:5n3)", inchikey = key),
                   list(name = "docosapentaenoate (n6 DPA; 22:5n6)", inchikey = key))
  expect_identical(r$rule_id, "R03c")
})

test_that("R03c also stops R06b", {
  r <- pair_result(
    list(name = "FA 22:5 n-3", lipid_species = "FA 22:5", lipid_level = "species"),
    list(name = "FA 22:5 n-6", lipid_species = "FA 22:5", lipid_level = "species"))
  expect_identical(r$rule_id, "R03c")
})

test_that("R03 stays exact when only one name states stereo, or none does", {
  key <- "QNAYBMKLOCPYGJ-REOHCLBHSA-N"
  r <- pair_result(list(name = "L-alanine", inchikey = key),
                   list(name = "alanine", inchikey = key))
  expect_identical(r$rule_id, "R03")
  key <- "DHMQDGOQFOQNFH-UHFFFAOYSA-N"
  r <- pair_result(list(name = "glycine", inchikey = key),
                   list(name = "Gly", inchikey = key))
  expect_identical(r$rule_id, "R03")
  key <- "BPYKTIZUTYGOLE-KDUUSRDASA-N"
  r <- pair_result(list(name = "bilirubin (Z,Z)", inchikey = key),
                   list(name = "bilirubin (Z,Z)*", inchikey = key))
  expect_identical(r$rule_id, "R03")
})

test_that("R99: keys that differ only in protonation are not called stereo variants", {
  r <- pair_result(list(name = "bilirubin (Z,Z)", inchikey = "BPYKTIZUTYGOLE-KDUUSRDASA-N"),
                   list(name = "bilirubin (Z,Z)", inchikey = "BPYKTIZUTYGOLE-KDUUSRDASA-M"))
  expect_identical(r$rule_id, "R99")
  expect_identical(r$relationship, "unresolved")
})

test_that("every pair gets exactly one known relationship and rule", {
  rows <- data.frame(
    row_id = c("a", "b", "c"),
    original_name = c("inosine", "Inosine", "Inosine-15N4 [iSTD]"),
    inchikey = c("UGQMRVRMYYASKQ-KQYNXXCUSA-N", "UGQMRVRMYYASKQ-KQYNXXCUSA-N", "inchikey"),
    lipid_species = NA_character_, lipid_molecular = NA_character_,
    lipid_level = NA_character_, group_id = "Inosine",
    stringsAsFactors = FALSE)
  out <- classify_pairs(row_evidence(rows), make_pairs(rows))
  expect_equal(nrow(out), 3)
  expect_false(anyNA(out$relationship))
  expect_true(all(out$relationship %in% metMatchAudit:::.relationships))
  expect_true(all(out$rule_id %in% metMatchAudit:::.rule_order))
})

test_that("missing evidence columns stop with a clear message", {
  rows <- row_evidence(data.frame(row_id = "a", original_name = "x"))
  expect_error(classify_pairs(rows, data.frame(group_id = "g", row_a = "a", row_b = "a")),
               "lipid_evidence")
})

test_that("R04 separates a group row from a different Species (the SM 36:1 review group)", {
  r <- pair_result(
    list(name = "sphingomyelin (d18:1/20:2, d18:2/20:1, d16:1/22:2)*",
         lipid_species = "SM 38:3;O2", lipid_level = "group"),
    list(name = "stearoyl sphingomyelin (d18:1/18:0)",
         lipid_species = "SM 36:1;O2", lipid_level = "structure"))
  expect_identical(r$rule_id, "R04")
  expect_identical(r$relationship, "distinct")
})

test_that("R05: a Species-level row is broader than a group row", {
  r <- pair_result(
    list(name = "sphingomyelin (d18:2/16:0, d18:1/16:1)*",
         lipid_species = "SM 34:2;O2", lipid_level = "group"),
    list(name = "SM 34:2;O2", lipid_species = "SM 34:2;O2", lipid_level = "species"))
  expect_identical(r$rule_id, "R05")
  expect_identical(r$relationship, "broader_narrower")
})

test_that("R05g: a group row against a structure stays unresolved", {
  # The structure may or may not be one of the listed members.
  r <- pair_result(
    list(name = "sphingomyelin (d18:2/16:0, d18:1/16:1)*",
         lipid_species = "SM 34:2;O2", lipid_level = "group"),
    list(name = "some SM with the same sum", lipid_species = "SM 34:2;O2",
         lipid_molecular = "SM 18:1;O2/16:1", lipid_level = "structure"))
  expect_identical(r$rule_id, "R05g")
  expect_identical(r$relationship, "unresolved")
})

test_that("R09b: partial stereo against full stereo is broader/narrower", {
  # Cholesterol, review label D. Without structures this was R11.
  r <- pair_result(
    list(name = "cholesterol", inchikey = "HVYWMOMLDIMFJA-LBHVWCRMSA-N",
         st_usable_inchi = inchi_cholesterol_leocc),
    list(name = "Cholesterol", inchikey = "HVYWMOMLDIMFJA-DPAQBDIFSA-N",
         st_usable_inchi = inchi_cholesterol_broad))
  expect_identical(r$rule_id, "R09b")
  expect_identical(r$relationship, "broader_narrower")
})

test_that("R09b does not fire on a stereo conflict; the pair goes on to R11", {
  r <- pair_result(
    list(name = "taurocholate", inchikey = "WBWWGRHZICKQGZ-HGVCTONRSA-N",
         st_usable_inchi = inchi_taurocholate_leocc),
    list(name = "Taurocholic acid", inchikey = "WBWWGRHZICKQGZ-HZAMXZRMSA-N",
         st_usable_inchi = inchi_taurocholate_broad))
  expect_identical(r$rule_id, "R11")
})

test_that("without trusted structures the pair stays at R11", {
  r <- pair_result(
    list(name = "cholesterol", inchikey = "HVYWMOMLDIMFJA-LBHVWCRMSA-N"),
    list(name = "Cholesterol", inchikey = "HVYWMOMLDIMFJA-DPAQBDIFSA-N",
         st_usable_inchi = inchi_cholesterol_broad))
  expect_identical(r$rule_id, "R11")
})
