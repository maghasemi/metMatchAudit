# RefMet check. PubChem is replaced by a fake that answers from a fixed
# list, so the tests run offline and show exactly which lookups happen.

refmet_text <- c(
  "name\tPhenyl acetate", "pubchem_cid\t31229", "inchi_key\tIPBVNPXQWQGGJP-UHFFFAOYSA-N",
  "formula\tC8H8O2", "refmet_id\tRM1", "",
  "name\tPhenylacetic acid", "pubchem_cid\t999", "inchi_key\tWLJVXDMOQOGPHL-UHFFFAOYSA-N",
  "formula\tC8H8O2", "refmet_id\tRM2", "",
  "name\tCholine", "pubchem_cid\t6209", "inchi_key\tSGMZJAMFUVOLNK-UHFFFAOYSA-M",
  "formula\tC5H14ClNO", "refmet_id\tRM3", "",
  "name\tFructose", "pubchem_cid\t2723872", "inchi_key\tLKDRXBCSQODPBY-VRPWFDPXSA-N",
  "formula\tC6H12O6", "refmet_id\tRM4", "",
  "name\tMG 18:2", "pubchem_cid\t5365676", "inchi_key\tWECGLUPZRHILCT-GSNKCQISSA-N",
  "formula\tC21H38O4", "refmet_id\tRM5"
)

refmet_db_for_tests <- function() {
  f <- tempfile(fileext = ".txt")
  writeLines(refmet_text, f)
  load_refmet(f)
}

# Parent CIDs and properties the fake PubChem knows about.
fake_pubchem <- function() {
  parents <- list(
    "inchikey/WLJVXDMOQOGPHL-UHFFFAOYSA-N" = "999",
    "inchikey/OEYIOHPDSNJKLS-UHFFFAOYSA-N" = "305",    # choline cation
    "cid/6209" = "305",                                 # choline chloride
    "inchikey/RFSUNEUAIZKAJO-ARQDHWQXSA-N" = "5984",    # fructofuranose
    "cid/2723872" = "2723872",
    "cid/31229" = "31229",
    "inchikey/IEPGNWMPIFDNSD-HZJYTTRNSA-N" = "5283469", # 2-linoleoylglycerol
    "inchikey/XXXXXXXXXXXXXX-UHFFFAOYSA-N" = "1"
  )
  props <- data.frame(
    CID = c(999, 305, 5984, 2723872, 31229, 5283469, 1),
    InChIKey = c("WLJVXDMOQOGPHL-UHFFFAOYSA-N", "OEYIOHPDSNJKLS-UHFFFAOYSA-N",
                 "RFSUNEUAIZKAJO-ARQDHWQXSA-N", "LKDRXBCSQODPBY-VRPWFDPXSA-N",
                 "IPBVNPXQWQGGJP-UHFFFAOYSA-N", "IEPGNWMPIFDNSD-HZJYTTRNSA-N",
                 "XXXXXXXXXXXXXX-UHFFFAOYSA-N"),
    MolecularFormula = c("C8H8O2", "C5H14NO+", "C6H12O6", "C6H12O6", "C8H8O2",
                         "C21H38O4", "C9H9O2"),
    stringsAsFactors = FALSE)
  calls <- character(0)
  fetch <- function(path, body = NULL) {
    calls <<- c(calls, path)
    m <- regmatches(path, regexec("^compound/(inchikey|cid)/([^/]+)/cids/JSON", path))[[1]]
    if (length(m) == 3) {
      p <- parents[[paste0(m[2], "/", m[3])]]
      if (is.null(p)) return(list(Fault = list(Code = "PUGREST.NotFound")))
      return(list(IdentifierList = list(CID = as.integer(p))))
    }
    if (startsWith(path, "compound/cid/property/")) {
      want <- strsplit(body$cid, ",")[[1]]
      return(list(PropertyTable = list(Properties = props[as.character(props$CID) %in% want, ])))
    }
    list(error = "unexpected path")
  }
  list(fetch = fetch, calls = function() calls)
}

check_rows <- function(name, group, inchikey = NA_character_, lipid_species = NA_character_,
                       fetch = fake_pubchem()$fetch) {
  rows <- data.frame(row_id = paste0("r", seq_along(name)), original_name = name,
                     group_id = group, inchikey = inchikey, stringsAsFactors = FALSE)
  ev <- row_evidence(rows)
  ev$lipid_species <- lipid_species
  refmet_check(ev, refmet_db_for_tests(), fetch = fetch)
}

test_that("load_refmet reads the text download and blanks invalid keys", {
  db <- refmet_db_for_tests()
  expect_s3_class(db, "refmet_db")
  expect_identical(nrow(db), 5L)
  expect_identical(db$formula[db$name == "Choline"], "C5H14ClNO")
  expect_true(nzchar(attr(db, "md5")))
  expect_error(load_refmet(tempfile()), "Cannot find")
})

test_that("a matching key needs no PubChem lookup", {
  pc <- fake_pubchem()
  out <- check_rows("phenylacetic acid", "Phenylacetic acid",
                    inchikey = "WLJVXDMOQOGPHL-UHFFFAOYSA-N", fetch = pc$fetch)
  expect_identical(out$rm_check, "agrees")
  expect_identical(out$rm_basis, "block1")
  expect_length(pc$calls(), 0)
})

test_that("phenylacetate assigned to phenyl acetate is an assignment conflict", {
  out <- check_rows("phenylacetate", "Phenyl acetate", inchikey = "WLJVXDMOQOGPHL-UHFFFAOYSA-N")
  expect_identical(out$rm_check, "conflict")
  expect_identical(out$rm_basis, "parent_block1")
})

test_that("a salt and its free form agree through their parent compound", {
  out <- check_rows("choline", "Choline", inchikey = "OEYIOHPDSNJKLS-UHFFFAOYSA-N")
  expect_identical(out$rm_check, "agrees")
  expect_identical(out$rm_basis, "parent_block1")
})

test_that("same name, different structure is an identifier conflict", {
  out <- check_rows("fructose", "Fructose", inchikey = "RFSUNEUAIZKAJO-ARQDHWQXSA-N")
  expect_identical(out$rm_check, "identifier_conflict")
})

test_that("an entry below structure level is compared by Species, not by its key", {
  out <- check_rows("2-linoleoylglycerol (18:2)", "MG 18:2",
                    inchikey = "IEPGNWMPIFDNSD-HZJYTTRNSA-N", lipid_species = "MG 18:2")
  expect_identical(out$rm_check, "agrees")
  expect_identical(out$rm_basis, "species")
})

test_that("an entry below structure level against a row with only a key compares formulas", {
  out <- check_rows("2-linoleoylglycerol (18:2)", "MG 18:2",
                    inchikey = "IEPGNWMPIFDNSD-HZJYTTRNSA-N")
  expect_identical(out$rm_check, "agrees")
  expect_identical(out$rm_basis, "formula")
  out <- check_rows("something else", "MG 18:2", inchikey = "XXXXXXXXXXXXXX-UHFFFAOYSA-N")
  expect_identical(out$rm_check, "conflict")
})

test_that("a missing PubChem answer leaves the row unchecked, not flagged", {
  out <- check_rows("phenylacetate", "Phenyl acetate", inchikey = "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N")
  expect_identical(out$rm_check, "unchecked")
})

test_that("group names, missing entries and rows without evidence are not compared", {
  out <- check_rows(c("mannitol/sorbitol", "x", "y", "z"),
                    c("Phenyl acetate", "Not in RefMet", "-", "Phenyl acetate"))
  expect_identical(out$rm_check, c("skipped", "entry_missing", "no_group", "not_compared"))
})

test_that("row_flags reports C06 and C07", {
  out <- check_rows(c("phenylacetate", "fructose"), c("Phenyl acetate", "Fructose"),
                    inchikey = c("WLJVXDMOQOGPHL-UHFFFAOYSA-N", "RFSUNEUAIZKAJO-ARQDHWQXSA-N"))
  out$lipid_conflict <- FALSE
  out$sl_species_abbrev <- NA_character_
  fl <- row_flags(out)
  expect_identical(fl$check, c("C06", "C07"))
  expect_identical(fl$flag, c("refmet_assignment_conflict", "refmet_identifier_conflict"))
})
