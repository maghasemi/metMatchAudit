# These need RDKit through reticulate and are skipped without it.

test_that("without RDKit use, the columns are added and left empty", {
  ev <- row_evidence(data.frame(original_name = "x", inchikey = NA))
  out <- structure_evidence(ev, use_rdkit = FALSE)
  expect_true(all(c("st_inchi", "st_usable_inchi") %in% names(out)))
  expect_true(is.na(out$st_usable_inchi))
})

test_that("a SMILES gives the key PubChem could not compute (LEOCC DHEA-S)", {
  skip_if_not(rdkit_available(), "RDKit not available")
  ev <- row_evidence(data.frame(
    original_name = "dehydroepiandrosterone sulfate (DHEA-S)",
    inchikey = "CZWCKYRVOZZJNM-DYIPZHMLSA-N",
    smiles = "[H][C@@]1(CC[C@]2(C)[C@]3(CCC2=O)[H])[C@@]4(C)C(C[C@H](OS(=O)(O)=O)CC4)=CC[C@]13[H]",
    stringsAsFactors = FALSE))
  out <- structure_evidence(ev)
  expect_identical(out$st_key, "CZWCKYRVOZZJNM-DYIPZHMLSA-N")
  expect_identical(out$st_source, "rdkit_smiles")
  expect_true(out$st_key_agrees)
  expect_false(is.na(out$st_usable_inchi))
})

test_that("a file InChI that disagrees with its key is not trusted (Broad xylose)", {
  skip_if_not(rdkit_available(), "RDKit not available")
  ev <- row_evidence(data.frame(
    original_name = "Xylose", inchikey = "SRBFZHDQGSBBOR-IOVATXLUSA-N",
    inchi = inchi_xylose_broad, stringsAsFactors = FALSE))
  out <- structure_evidence(ev)
  expect_identical(out$st_key, "PYMYPHUHKUWMLA-VPENINKCSA-N")
  expect_true(out$st_block1_differs)
  expect_true(is.na(out$st_usable_inchi))
})
