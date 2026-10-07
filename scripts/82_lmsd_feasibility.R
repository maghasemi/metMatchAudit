# Feasibility for a second lipid source: LIPID MAPS Structure Database
# (LMSD), looked up by InChIKey the same way SwissLipids is.
#
# Counts only. No expert label is read and no package code changes.
# Criteria L1 to L3 are in notes/decisions.md (2026-10-07, LMSD plan).
#
# Needs: data-raw/lipidmaps/structures.sdf (LMSD standard SDF, 2026-03-06),
# results/biocrates_feasibility/table_mapped.csv from scripts/81, and the
# v0.4.2 audits of the development data and of COMETS.
# Output: data-raw/lipidmaps/lmsd_slim.rds, results/lmsd_feasibility/
# Run from the project root: Rscript scripts/82_lmsd_feasibility.R

pkgload::load_all(".", quiet = TRUE)

sdf_path     <- "data-raw/lipidmaps/structures.sdf"
expected_md5 <- "134c2b5c6d76e9883d46339a6289eefd"
slim_path    <- "data-raw/lipidmaps/lmsd_slim.rds"
out_dir      <- "results/lmsd_feasibility"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (unname(tools::md5sum(sdf_path)) != expected_md5) {
  stop("md5 of ", sdf_path, " differs from the file this plan was written for.")
}


# ---- Read the SDF once, keep the fields we need -------------------------------
# Each field is a "> <NAME>" line followed by its value. Records end at "$$$$".

if (!file.exists(slim_path)) {
  x <- readLines(sdf_path, warn = FALSE, encoding = "UTF-8")
  rec <- cumsum(c(0L, head(x == "$$$$", -1))) + 1L
  n <- max(rec)
  fields <- c("LM_ID", "ABBREVIATION", "INCHI_KEY", "MAIN_CLASS", "SUB_CLASS")
  lmsd <- data.frame(matrix(NA_character_, n, length(fields), dimnames = list(NULL, tolower(fields))),
                     stringsAsFactors = FALSE)
  for (f in fields) {
    i <- which(x == paste0("> <", f, ">"))
    lmsd[[tolower(f)]][rec[i]] <- x[i + 1]
  }
  rm(x)
  lmsd <- lmsd[!is.na(lmsd$lm_id), ]
  lmsd$inchi_key[!is_valid_inchikey(lmsd$inchi_key)] <- NA_character_
  saveRDS(lmsd, slim_path)
}
lmsd <- readRDS(slim_path)
message("LMSD records: ", nrow(lmsd), ", with a valid InChIKey: ", sum(!is.na(lmsd$inchi_key)),
        ", with an abbreviation: ", sum(!is.na(lmsd$abbreviation)))


# ---- Lookup, fixed in the plan --------------------------------------------------
# Full key, then the first 25 characters, then the first 14. Stop at the
# first tier with a hit. The row gets a Species only if every hit's
# abbreviation parses and they all give the same Species.

k <- lmsd[!is.na(lmsd$inchi_key) & !is.na(lmsd$abbreviation), ]
k$k25 <- substr(k$inchi_key, 1, 25)
k$k14 <- substr(k$inchi_key, 1, 14)

by_full <- split(k$abbreviation, k$inchi_key)
by_k25  <- split(k$abbreviation, k$k25)
by_k14  <- split(k$abbreviation, k$k14)

lookup <- function(keys) {
  hits <- lapply(keys, function(key) {
    if (is.na(key)) return(list(tier = NA_character_, abbr = character(0)))
    tiers <- list(full = by_full[[key]], k25 = by_k25[[substr(key, 1, 25)]],
                  k14 = by_k14[[substr(key, 1, 14)]])
    for (tier in names(tiers)) {
      a <- unique(tiers[[tier]])
      if (length(a) > 0) return(list(tier = tier, abbr = a))
    }
    list(tier = NA_character_, abbr = character(0))
  })
  all_abbr <- unique(unlist(lapply(hits, `[[`, "abbr")))
  lp <- parse_lipid_names(all_abbr)
  sp <- stats::setNames(ifelse(lp$parsed, lp$species, NA_character_), all_abbr)
  data.frame(
    lmsd_tier = vapply(hits, `[[`, character(1), "tier"),
    lmsd_hits = vapply(hits, function(h) length(h$abbr), integer(1)),
    lmsd_species = vapply(hits, function(h) {
      s <- unique(sp[h$abbr])
      if (length(s) == 1 && !is.na(s)) s else NA_character_
    }, character(1)),
    stringsAsFactors = FALSE)
}

check <- function(audit, label) {
  r <- audit$rows
  eligible <- r$ik_valid & !r$name_group & !r$name_alternatives & !r$name_internal_standard
  res <- lookup(ifelse(eligible, r$inchikey, NA_character_))
  r <- cbind(r, res)

  had <- !is.na(r$lipid_species)
  gain <- eligible & !had & !is.na(r$lmsd_species)
  both <- had & !is.na(r$lmsd_species)
  agree <- both & (r$lipid_species == r$lmsd_species |
                     .generic_headgroup(r$lipid_species) == .generic_headgroup(r$lmsd_species))
  message("\n==== ", label)
  message("Rows with a key that can be looked up: ", sum(eligible))
  message("L1, rows that get a Species from LMSD and have none today: ", sum(gain))
  print(table(r$file[gain]))
  message("Hit tier for those rows:")
  print(table(r$lmsd_tier[gain]))
  message("L3, rows with a Species from both sources: ", sum(both), ", agree: ", sum(agree),
          if (sum(both) > 0) sprintf(" (%.1f%%)", 100 * sum(agree) / sum(both)) else "")
  dis <- r[both & !agree, c("file", "original_name", "lipid_species", "lmsd_species", "lmsd_tier")]
  if (nrow(dis) > 0) { message("Disagreements:"); print(dis, row.names = FALSE, right = FALSE) }
  write.csv(r[gain, c("row_id", "file", "original_name", "lmsd_species", "lmsd_tier", "lmsd_hits")],
            file.path(out_dir, paste0(label, "_gained.csv")), row.names = FALSE)
  r
}

dev    <- check(readRDS("results/audit_v0.1_sealed.rds"), "development")
comets_audit <- readRDS("results/comets_v0.4.2/audit.rds")
comets <- check(comets_audit, "comets")


# ---- L2: the F3 count again, with both sources ------------------------------------

tab <- read.csv("results/biocrates_feasibility/table_mapped.csv", colClasses = "character",
                na.strings = c("", "NA"))
bio <- comets$file == "inputs_1539353822"
abbr <- regmatches(comets$original_name,
                   regexec("\\(\\s*((?:[^()]|\\([^()]*\\))*?)\\s*\\)\\s*$", comets$original_name, perl = TRUE))
key <- vapply(abbr, function(v) if (length(v) == 2) gsub("\\s+", "", v[2]) else NA_character_, character(1))
bio_sp <- ifelse(bio, tab$species[match(key, tab$key)], NA_character_)

sp_any <- ifelse(!is.na(comets$lipid_species), comets$lipid_species, comets$lmsd_species)
id <- comets$row_id
bs <- stats::setNames(bio_sp, id); os <- stats::setNames(sp_any, id); fl <- stats::setNames(comets$file, id)
p <- comets_audit$pairs
p <- p[fl[p$row_a] != fl[p$row_b] & (!is.na(bs[p$row_a]) | !is.na(bs[p$row_b])), ]
other <- ifelse(!is.na(bs[p$row_a]), p$row_b, p$row_a)
message("\nL2, cross-file pairs with a mapped Biocrates row: ", nrow(p),
        "; other row has a Species from SwissLipids or LMSD: ", sum(!is.na(os[other])))
