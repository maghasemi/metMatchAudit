# PubChem is replaced by a fake that answers from fixed tables, so these
# tests run offline. Keys and CIDs are real: choline chloride (CID 6209)
# has parent choline (CID 305); inosine is CID 6021.
#
# `calls$fail` is a pattern: requests whose path matches it get the busy
# answer PubChem sends under load. Set it to NULL to let them through.
busy <- list(Fault = list(Code = "PUGREST.ServerBusy", Message = "Too many requests or server too busy"))
fake_pubchem <- function(calls = new.env(), fail = NULL) {
  calls$n <- 0
  calls$fail <- fail
  parents <- c("6209" = "305", "6021" = "6021", "111" = "111", "222" = "222")
  keys <- c("305"  = "OEYIOHPDSNJKLS-UHFFFAOYSA-N",
            "6021" = "UGQMRVRMYYASKQ-KQYNXXCUSA-N",
            "111"  = "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
            "222"  = "BBBBBBBBBBBBBB-UHFFFAOYSA-N")
  hmdb <- list(HMDB0000195 = "6021", HMDB0000097 = c("6209", "305"),
               HMDB0009999 = c("111", "222"))
  fetch <- function(path, body = NULL) {
    calls$n <- calls$n + 1
    if (!is.null(calls$fail) && grepl(calls$fail, path)) return(busy)
    if (grepl("cids_type=parent", path)) {
      cid <- sub("compound/cid/([0-9]+)/.*", "\\1", path)
      return(list(IdentifierList = list(CID = unname(parents[cid]))))
    }
    if (grepl("^compound/cid/property", path)) {
      cids <- strsplit(body$cid, ",")[[1]]
      cids <- cids[cids %in% names(keys)]
      return(list(PropertyTable = list(Properties = data.frame(CID = cids, InChIKey = unname(keys[cids])))))
    }
    id <- sub(".*RegistryID/(HMDB[0-9]+)/.*", "\\1", path)
    if (is.null(hmdb[[id]])) return(list(Fault = list(Code = "PUGREST.NotFound", Message = "No CID found")))
    list(IdentifierList = list(CID = hmdb[[id]]))
  }
  list(fetch = fetch, calls = calls)
}

test_that("identifier cells are normalised", {
  expect_identical(normalize_pubchem_id(c("pubchem:96", "12594", "pubchem:439685,7405", NA, "abc")),
                   c("96", "12594", NA, NA, NA))
  expect_identical(normalize_hmdb_id(c("HMDB00641", "hmdb:HMDB14493", "hmdb:HMDB0240679,HMDB0240680", NA)),
                   c("HMDB0000641", "HMDB0014493", NA, NA))
})

rows_for_resolution <- function() {
  data.frame(
    original_name = c("choline", "inosine", "glutamine", "Citric and or Isocitric acid",
                      "some compound", "two ids", "unknown"),
    inchikey = c(NA, NA, "ZDXPYRJPNDTMRX-VKHMYHEASA-N", NA, NA, NA, NA),
    pubchem  = c("6209", NA, NA, "311", NA, "pubchem:6209,6021", NA),
    hmdb     = c(NA, "HMDB00195", NA, NA, "HMDB0009999", NA, NA),
    stringsAsFactors = FALSE)
}

test_that("a CID gives the key of its parent compound", {
  fk <- fake_pubchem()
  out <- resolve_identifiers(rows_for_resolution(), fetch = fk$fetch)
  # Choline chloride resolves to choline itself.
  expect_identical(out$inchikey[1], "OEYIOHPDSNJKLS-UHFFFAOYSA-N")
  expect_identical(out$inchikey_source[1], "pubchem_cid")
})

test_that("an HMDB ID in the old five-digit form resolves", {
  out <- resolve_identifiers(rows_for_resolution(), fetch = fake_pubchem()$fetch)
  expect_identical(out$inchikey[2], "UGQMRVRMYYASKQ-KQYNXXCUSA-N")
  expect_identical(out$inchikey_source[2], "hmdb_via_pubchem")
})

test_that("a key from the file is kept and marked as such", {
  out <- resolve_identifiers(rows_for_resolution(), fetch = fake_pubchem()$fetch)
  expect_identical(out$inchikey[3], "ZDXPYRJPNDTMRX-VKHMYHEASA-N")
  expect_identical(out$inchikey_source[3], "file")
})

test_that("group names, lists, and HMDB IDs with several keys stay unresolved", {
  out <- resolve_identifiers(rows_for_resolution(), fetch = fake_pubchem()$fetch)
  expect_true(all(is.na(out$inchikey[c(4, 5, 6, 7)])))
  expect_true(all(is.na(out$inchikey_source[c(4, 5, 6, 7)])))
})

test_that("an invalid key in the file is kept in inchikey_file and replaced", {
  r <- data.frame(original_name = "Inosine-like", inchikey = "inchikey", pubchem = "6021")
  out <- resolve_identifiers(r, fetch = fake_pubchem()$fetch)
  expect_identical(out$inchikey_file, "inchikey")
  expect_identical(out$inchikey, "UGQMRVRMYYASKQ-KQYNXXCUSA-N")
})

test_that("answers are cached on disk and reused", {
  f <- tempfile(fileext = ".rds")
  fk <- fake_pubchem()
  resolve_identifiers(rows_for_resolution(), cache_file = f, fetch = fk$fetch)
  first <- fk$calls$n
  expect_gt(first, 0)
  resolve_identifiers(rows_for_resolution(), cache_file = f, fetch = fk$fetch)
  expect_identical(fk$calls$n, first)
})

test_that("rows that need nothing never call PubChem", {
  fk <- fake_pubchem()
  r <- data.frame(original_name = "glutamine", inchikey = "ZDXPYRJPNDTMRX-VKHMYHEASA-N")
  out <- resolve_identifiers(r, fetch = fk$fetch)
  expect_identical(fk$calls$n, 0)
  expect_identical(out$inchikey_source, "file")
})

test_that("C01 looks at the value from the file, not the resolved key", {
  r <- data.frame(row_id = "x", original_name = "Inosine-like", inchikey = "inchikey", pubchem = "6021")
  ev <- row_evidence(resolve_identifiers(r, fetch = fake_pubchem()$fetch))
  ev$lipid_conflict <- FALSE; ev$sl_species_abbrev <- NA
  f <- row_flags(ev)
  expect_identical(f$check, "C01")
})

test_that("pair evidence names a key that came from PubChem", {
  rows <- data.frame(row_id = c("a", "b"), group_id = "Inosine",
                     original_name = c("inosine", "Inosine"),
                     inchikey = c("UGQMRVRMYYASKQ-KQYNXXCUSA-N", NA),
                     pubchem = c(NA, "6021"), stringsAsFactors = FALSE)
  ev <- row_evidence(resolve_identifiers(rows, fetch = fake_pubchem()$fetch))
  ev$lipid_species <- NA; ev$lipid_molecular <- NA; ev$lipid_level <- NA
  p <- classify_pairs(ev, make_pairs(ev))
  expect_identical(p$rule_id, "R03")
  expect_match(p$evidence, "second key from pubchem_cid")
})

test_that("Species-level lipid names are not resolved, structure-level ones are", {
  # "TG 42:0" names a sum composition, so an identifier on that row can
  # only point to one representative structure.
  fk <- fake_pubchem()
  r <- data.frame(original_name = c("TG 42:0", "SM 36:1;O2", "inosine"),
                  inchikey = NA, pubchem = c("111", "222", "6021"),
                  stringsAsFactors = FALSE)
  out <- resolve_identifiers(r, fetch = fk$fetch)
  expect_true(all(is.na(out$inchikey[1:2])))
  expect_identical(out$inchikey[3], "UGQMRVRMYYASKQ-KQYNXXCUSA-N")
})

test_that("internal standards are not resolved", {
  # Broad "Inosine-15N4 [iSTD]" carries the HMDB ID of plain inosine.
  r <- data.frame(original_name = "Inosine-15N4 [iSTD]", inchikey = "inchikey",
                  hmdb = "HMDB0000195", stringsAsFactors = FALSE)
  out <- resolve_identifiers(r, fetch = fake_pubchem()$fetch)
  expect_identical(out$inchikey, "inchikey")
  expect_true(is.na(out$inchikey_source))
})

# ---- Failed requests (0.9.1) --------------------------------------------------

choline_row <- function() data.frame(original_name = "choline", inchikey = NA, pubchem = "6209")

test_that("a failed parent request leaves the row without a key and is counted", {
  # Before 0.9.1 the CID itself was used, so the row got the key of
  # choline chloride's own record as if it were an answer.
  fk <- fake_pubchem(fail = "cids_type=parent")
  expect_message(out <- resolve_identifiers(choline_row(), fetch = fk$fetch), "1 PubChem request failed")
  expect_true(is.na(out$inchikey))
  expect_true(is.na(out$inchikey_source))
  expect_identical(attr(out, "pubchem_failed"), 1)
})

test_that("a failed key request is counted and gives no key", {
  fk <- fake_pubchem(fail = "^compound/cid/property")
  expect_message(out <- resolve_identifiers(choline_row(), fetch = fk$fetch), "PubChem request failed")
  expect_true(is.na(out$inchikey))
  expect_identical(attr(out, "pubchem_failed"), 1)
})

test_that("a busy answer is not cached, so the next run asks again", {
  f <- tempfile(fileext = ".rds")
  fk <- fake_pubchem(fail = "cids_type=parent")
  expect_message(resolve_identifiers(choline_row(), cache_file = f, fetch = fk$fetch), "failed")
  fk$calls$fail <- NULL
  expect_silent(out <- resolve_identifiers(choline_row(), cache_file = f, fetch = fk$fetch))
  expect_identical(out$inchikey, "OEYIOHPDSNJKLS-UHFFFAOYSA-N")
  expect_identical(attr(out, "pubchem_failed"), 0)
})

test_that("a busy answer already in a cache is asked again", {
  # Versions before 0.9.1 kept such answers.
  f <- tempfile(fileext = ".rds")
  env <- new.env()
  env[[paste("compound/cid/6209/cids/JSON?cids_type=parent", "")]] <- busy
  saveRDS(env, f)
  fk <- fake_pubchem()
  out <- resolve_identifiers(choline_row(), cache_file = f, fetch = fk$fetch)
  expect_identical(out$inchikey, "OEYIOHPDSNJKLS-UHFFFAOYSA-N")
})

test_that("not found is a real answer: cached, and not counted as failed", {
  f <- tempfile(fileext = ".rds")
  r <- data.frame(original_name = "unknown thing", inchikey = NA, hmdb = "HMDB0001234")
  fk <- fake_pubchem()
  expect_silent(out <- resolve_identifiers(r, cache_file = f, fetch = fk$fetch))
  expect_identical(attr(out, "pubchem_failed"), 0)
  first <- fk$calls$n
  resolve_identifiers(r, cache_file = f, fetch = fk$fetch)
  expect_identical(fk$calls$n, first)
})

test_that("an HMDB row whose CIDs could not all be asked stays open", {
  # HMDB0000097 lists CIDs 6209 and 305, both leading to choline. If the
  # parent of 6209 cannot be asked, the one answer left is not enough.
  r <- data.frame(original_name = "choline", inchikey = NA, hmdb = "HMDB0000097")
  fk <- fake_pubchem(fail = "compound/cid/6209/")
  expect_message(out <- resolve_identifiers(r, fetch = fk$fetch), "failed")
  expect_true(is.na(out$inchikey))
  fk$calls$fail <- NULL
  out <- resolve_identifiers(r, fetch = fk$fetch)
  expect_identical(out$inchikey, "OEYIOHPDSNJKLS-UHFFFAOYSA-N")
})

test_that("the message says how to keep answers when there is no cache file", {
  fk <- fake_pubchem(fail = "cids_type=parent")
  expect_message(resolve_identifiers(choline_row(), fetch = fk$fetch), "cache_file")
})
