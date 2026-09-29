# Lipid feasibility test. The criteria were written into
# notes/decisions.md before this script was run for the first time.
#
# Two questions:
#   1. Can each Metabolon (LEOCC) lipid row be placed in the SwissLipids
#      hierarchy through its InChIKey, so we know its Species-level name?
#   2. Can each Broad lipid name be parsed by rgoslin to Species level
#      or finer?
# Plus a check against the manual review: for lipid rows in groups
# already labeled in notes/review.csv, print what SwissLipids and
# rgoslin say, so it can be compared with the manual label by eye.
#
# Needs: scripts/10_swisslipids_slim.R run first, and the MetLinkR data
# from scripts/00_get_data.sh.
#
# Output folder: results/lipid_feasibility/
# Run from the project root: Rscript scripts/11_lipid_feasibility.R

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
})

if (!requireNamespace("rgoslin", quietly = TRUE)) {
  stop("rgoslin is not installed. Run renv::install(\"bioc::rgoslin\").")
}

out_dir <- "results/lipid_feasibility"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- Scope, fixed in notes/decisions.md before the first run ---------

# LEOCC: rows with SUPER_PATHWAY "Lipid". Complex lipids are these
# SUB_PATHWAY values.
leocc_complex <- c("Phosphatidylcholine (PC)", "Phosphatidylethanolamine (PE)",
                   "Phosphatidylinositol (PI)", "Phosphatidylserine (PS)",
                   "Lysophospholipid", "Plasmalogen", "Lysoplasmalogen",
                   "Diacylglycerol", "Monoacylglycerol",
                   "Sphingomyelins", "Dihydrosphingomyelins",
                   "Ceramides", "Dihydroceramides",
                   "Hexosylceramides (HCER)", "Lactosylceramides (LCER)",
                   "Ceramide PEs")

# Broad: rows whose superClass is one of the eight LIPID MAPS
# categories. Complex lipids are these mainClass values.
lipidmaps_categories <- c("Fatty Acyls", "Glycerolipids", "Glycerophospholipids",
                          "Sphingolipids", "Sterol Lipids", "Prenol Lipids",
                          "Saccharolipids", "Polyketides")
broad_complex <- c("Glycerophosphocholines", "Glycerophosphoethanolamines",
                   "Monoradylglycerols", "Diradylglycerols", "Triradylglycerols",
                   "Sphingomyelins", "Ceramides")


# ---- Helpers ----------------------------------------------------------

key_ok <- function(x) !is.na(x) & grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", x)

# A group name reports several possible compounds in one row, for
# example "sphingomyelin (d18:1/20:2, d18:2/20:1)*" or
# "(2 or 3)-decenoate". Its InChIKey is one member picked from the
# group, so we do not look it up.
is_group_name <- function(x) {
  grepl(" or |and/or|and or| \\+ ", x, ignore.case = TRUE) |
    grepl("\\([^()]*:[^()]*,[^()]*\\)", x)
}

# rgoslin fills some not-applicable fields with the text "NA" rather
# than a real missing value. Turn those into real NA.
fix_na <- function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "")] <- NA_character_
  x
}

parse_fields <- c("Normalized.Name", "Grammar", "Level", "Species.Name",
                  "Molecular.Species.Name", "Sn.Position.Name", "Sum.Formula")

# Parse one name at a time. In our first trial, one unparseable name in
# a batch produced an odd error message, so each call stands alone.
parse_one <- function(x) {
  empty <- as.list(setNames(rep(NA_character_, length(parse_fields)), parse_fields))
  if (is.na(x) || !nzchar(x)) return(empty)
  r <- NULL
  invisible(capture.output(
    r <- tryCatch(suppressMessages(rgoslin::parseLipidNames(x)),
                  error = function(e) NULL)
  ))
  if (is.null(r) || nrow(r) == 0) return(empty)
  out <- lapply(parse_fields, function(f) {
    if (f %in% names(r)) fix_na(r[[f]][1]) else NA_character_
  })
  setNames(out, parse_fields)
}

parse_many <- function(x) {
  u <- unique(x[!is.na(x) & nzchar(x)])
  parsed <- rbindlist(lapply(u, function(v) as.data.table(parse_one(v))))
  parsed[, input := u]
  parsed
}

# Species level or finer, in rgoslin's terms.
finer_than_class <- function(level) {
  !is.na(level) & !(level %in% c("CATEGORY", "CLASS"))
}

pct <- function(a, b) if (b == 0) NA_real_ else round(100 * a / b, 1)


# ---- Load data ----------------------------------------------------------

sl <- readRDS("data-raw/swisslipids/lipids_slim.rds")
setDT(sl)

needed_levels <- c("Species", "Molecular subspecies",
                   "Structural subspecies", "Isomeric subspecies")
missing_levels <- setdiff(needed_levels, unique(sl$level))
if (length(missing_levels) > 0) {
  print(sort(unique(sl$level)))
  stop("Level names differ from what this script expects: ",
       paste(missing_levels, collapse = ", "))
}

base <- "data-raw/MetLinkR/inst/extdata/Example_output"
if (!dir.exists(base)) stop("Cannot find ", base, ". Run scripts/00_get_data.sh first.")

read_chr <- function(path) {
  read_csv(path, col_types = cols(.default = col_character()),
           locale = locale(encoding = "latin1"), na = c("", "NA"),
           show_col_types = FALSE)
}
leocc <- read_chr(file.path(base, "LEOCC_Metabolon_Annotations_metLinkR.csv"))
broad <- read_chr(file.path(base, "Broad_2022Aug_annotations_metLinkR.csv"))
stopifnot(nrow(leocc) == 1548, nrow(broad) == 511)

# row_id follows the same pattern as scripts/01_build_rows.R, so rows
# can be joined to review.csv.
rows <- rbind(
  data.table(row_id       = paste0("LEOCC_Metabolon:", seq_len(nrow(leocc))),
             file         = "LEOCC_Metabolon",
             name         = leocc$CHEMICAL_NAME,
             inchikey     = trimws(leocc$INCHIKEY),
             source_class = leocc$SUB_PATHWAY,
             is_lipid     = leocc$SUPER_PATHWAY %in% "Lipid",
             is_complex   = leocc$SUB_PATHWAY %in% leocc_complex &
                            leocc$SUPER_PATHWAY %in% "Lipid"),
  data.table(row_id       = paste0("Broad_2022Aug:", seq_len(nrow(broad))),
             file         = "Broad_2022Aug",
             name         = broad$Metabolite,
             inchikey     = trimws(broad$inChIKey),
             source_class = broad$mainClass,
             is_lipid     = broad$superClass %in% lipidmaps_categories,
             is_complex   = broad$mainClass %in% broad_complex &
                            broad$superClass %in% lipidmaps_categories)
)
rows <- rows[is_lipid == TRUE]
rows[, group_name := is_group_name(name)]
rows[, key_valid := key_ok(inchikey)]


# ---- SwissLipids lookup by InChIKey, in three tiers --------------------
# Tier 1: full key. Tier 2: first 25 characters, which ignores the
# protonation letter (SwissLipids stores structures as charged at
# pH 7.3). Tier 3: first 14 characters, the connectivity block. A row
# stops at the first tier that finds anything.

sl_keys <- sl[!is.na(inchikey), .(slm_id, full = inchikey, key25, key14)]

# Group-name rows are looked up too, but only to measure something:
# Metabolon reports many sphingomyelins as a group of molecular species
# that share one Species-level sum. Their results are kept apart and do
# not count toward the criteria.
query <- rows[key_valid == TRUE,
              .(row_id, full = inchikey,
                key25 = substr(inchikey, 1, 25),
                key14 = substr(inchikey, 1, 14))]

empty_hits <- function() {
  data.table(row_id = character(), slm_id = character(), tier = integer())
}

match_tier <- function(q, col, tier_value) {
  if (nrow(q) == 0) return(empty_hits())
  m <- merge(q[, .(row_id, k = get(col))], sl_keys[, .(slm_id, k = get(col))],
             by = "k", allow.cartesian = TRUE)
  # Check for zero rows before adding the tier column. Otherwise a
  # length-one tier value next to empty columns becomes a fake NA row.
  if (nrow(m) == 0) return(empty_hits())
  m[, .(row_id, slm_id, tier = tier_value)]
}

hits1 <- match_tier(query, "full", 1L)
left  <- query[!row_id %in% hits1$row_id]
hits2 <- match_tier(left, "key25", 2L)
left  <- left[!row_id %in% hits2$row_id]
hits3 <- match_tier(left, "key14", 3L)
hits  <- unique(rbind(hits1, hits2, hits3))


# ---- Walk up the hierarchy ----------------------------------------------
# An environment gives fast lookups by ID over 780k rows.

parent_env <- list2env(as.list(setNames(sl$parent, sl$slm_id)))
level_env  <- list2env(as.list(setNames(sl$level,  sl$slm_id)))
abbrev_env <- list2env(as.list(setNames(sl$abbrev, sl$slm_id)))

get_env <- function(ids, env) {
  if (length(ids) == 0) return(character(0))
  unlist(mget(ids, envir = env, ifnotfound = list(NA_character_)), use.names = FALSE)
}

split_ids <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  ids <- trimws(unlist(strsplit(x, "[|,]")))
  ids[grepl("^SLM:", ids)]
}

# All ancestors of one ID. Stops at 15 steps and never revisits an ID,
# so a loop in the parent links cannot hang the script.
ancestors <- function(id) {
  seen <- character(0)
  frontier <- id
  for (step in 1:15) {
    nxt <- setdiff(split_ids(get_env(frontier, parent_env)), c(seen, id))
    if (length(nxt) == 0) break
    seen <- c(seen, nxt)
    frontier <- nxt
  }
  seen
}

collapse <- function(x) {
  x <- sort(unique(x[!is.na(x)]))
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}

hit_info <- rbindlist(lapply(unique(hits$slm_id), function(id) {
  # Include the hit itself, in case it already sits at one of the
  # levels we want (for example a Species entry).
  anc <- c(id, ancestors(id))
  anc_level <- get_env(anc, level_env)
  data.table(slm_id           = id,
             hit_level        = get_env(id, level_env),
             hit_abbrev       = get_env(id, abbrev_env),
             species_abbrev   = collapse(get_env(anc[anc_level == "Species"], abbrev_env)),
             molecular_abbrev = collapse(get_env(anc[anc_level == "Molecular subspecies"], abbrev_env)),
             structural_abbrev = collapse(get_env(anc[anc_level == "Structural subspecies"], abbrev_env)))
}))

if (nrow(hit_info) == 0) {
  hit_info <- data.table(slm_id = character(), hit_level = character(),
                         hit_abbrev = character(), species_abbrev = character(),
                         molecular_abbrev = character(), structural_abbrev = character())
}

hits <- merge(hits, hit_info, by = "slm_id", all.x = TRUE)

# One line per row. A row can hit several SwissLipids entries, mostly
# in tier 3. Values are collapsed so nothing is silently dropped.
per_row <- hits[, .(tier              = min(tier),
                    n_hits            = uniqueN(slm_id),
                    hit_levels        = collapse(hit_level),
                    hit_abbrevs       = collapse(hit_abbrev),
                    species_abbrev    = collapse(unlist(strsplit(species_abbrev, "; "))),
                    molecular_abbrev  = collapse(unlist(strsplit(molecular_abbrev, "; "))),
                    structural_abbrev = collapse(unlist(strsplit(structural_abbrev, "; ")))),
                by = row_id]
per_row[, n_species := ifelse(is.na(species_abbrev), 0L,
                              lengths(strsplit(species_abbrev, "; ")))]

rows <- merge(rows, per_row, by = "row_id", all.x = TRUE)
rows[is.na(n_hits), n_hits := 0L]
rows[is.na(n_species), n_species := 0L]

# ---- Connectivity fallback (added after the first run) ------------------
# In the first run, some rows matched a SwissLipids entry that has no
# Level and no parent (three 1-monoacylglycerols), or an entry whose
# parent chain never reaches a Species (two d18:2 ceramides). For those
# rows only, we look at every SwissLipids entry with the same first 14
# InChIKey characters. Same connectivity means same constitution, so
# those entries must share one Species. We accept the fallback only if
# they all point to exactly one Species.
#
# The fallback fills new columns and leaves the original ones alone,
# so the pre-registered criteria can still be read from the original
# columns. See notes/decisions.md, 2026-09-29.

need_fallback <- rows[n_hits > 0 & n_species == 0 & key_valid == TRUE,
                      .(row_id, key14 = substr(inchikey, 1, 14))]

rows[, `:=`(species_fallback = NA_character_, n_species_fallback = 0L)]

if (nrow(need_fallback) > 0) {
  fb_hits <- merge(need_fallback, sl_keys[, .(slm_id, key14)],
                   by = "key14", allow.cartesian = TRUE)
  if (nrow(fb_hits) > 0) {
    fb_species <- fb_hits[, {
      sp <- unlist(lapply(slm_id, function(id) {
        anc <- c(id, ancestors(id))
        get_env(anc[get_env(anc, level_env) == "Species"], abbrev_env)
      }))
      sp <- sort(unique(sp[!is.na(sp)]))
      list(species_fallback   = if (length(sp) == 0) NA_character_ else paste(sp, collapse = "; "),
        n_species_fallback = length(sp))
    }, by = row_id]
    rows[fb_species, on = "row_id",
         `:=`(species_fallback = i.species_fallback,
              n_species_fallback = i.n_species_fallback)]
  }
}


# ---- rgoslin ---------------------------------------------------------------
# Three things get parsed:
#   the row's own name (works for shorthand names like Broad's),
#   the SwissLipids Species abbreviation (puts Metabolon rows into the
#     same notation as Broad),
#   the SwissLipids abbreviation of the hit itself (a second route to the
#     Species name, used as a cross-check against the hierarchy walk).

name_parsed <- parse_many(rows$name)
setnames(name_parsed, parse_fields, paste0("name_", parse_fields))
rows <- merge(rows, name_parsed, by.x = "name", by.y = "input", all.x = TRUE, sort = FALSE)

single_species <- rows[n_species == 1, species_abbrev]
species_parsed <- parse_many(single_species)
if (nrow(species_parsed) > 0) {
  rows <- merge(rows,
                species_parsed[, .(input, species_from_hierarchy = Species.Name)],
                by.x = "species_abbrev", by.y = "input", all.x = TRUE, sort = FALSE)
} else {
  rows[, species_from_hierarchy := NA_character_]
}

single_hit <- rows[n_hits == 1, hit_abbrevs]
hit_parsed <- parse_many(single_hit)
if (nrow(hit_parsed) > 0) {
  rows <- merge(rows,
                hit_parsed[, .(input, species_from_hit = Species.Name,
                               molecular_from_hit = Molecular.Species.Name)],
                by.x = "hit_abbrevs", by.y = "input", all.x = TRUE, sort = FALSE)
} else {
  rows[, `:=`(species_from_hit = NA_character_, molecular_from_hit = NA_character_)]
}

fb_single <- rows[n_species_fallback == 1, species_fallback]
fb_parsed <- parse_many(fb_single)
if (nrow(fb_parsed) > 0) {
  rows <- merge(rows,
                fb_parsed[, .(input, species_from_fallback = Species.Name)],
                by.x = "species_fallback", by.y = "input", all.x = TRUE, sort = FALSE)
} else {
  rows[, species_from_fallback := NA_character_]
}

setorder(rows, file, row_id)
fwrite(rows, file.path(out_dir, "lipid_rows.csv"))


# ---- Coverage tables ----------------------------------------------------

leocc_cov <- rows[file == "LEOCC_Metabolon",
                  .(rows          = .N,
                    group_names   = sum(group_name),
                    invalid_key   = sum(!key_valid),
                    found_tier1   = sum(!group_name & tier %in% 1L),
                    found_tier2   = sum(!group_name & tier %in% 2L),
                    found_tier3   = sum(!group_name & tier %in% 3L),
                    not_found     = sum(!group_name & key_valid & n_hits == 0),
                    one_species   = sum(!group_name & n_species == 1),
                    many_species  = sum(!group_name & n_species > 1),
                    species_parsed = sum(!group_name & !is.na(species_from_hierarchy)),
                    group_one_species = sum(group_name & n_species == 1)),
                  by = .(complex = is_complex, source_class)][order(-complex, -rows)]

broad_cov <- rows[file == "Broad_2022Aug",
                  .(rows            = .N,
                    parsed_species_or_finer = sum(finer_than_class(name_Level)),
                    found_in_swisslipids    = sum(n_hits > 0)),
                  by = .(complex = is_complex, source_class)][order(-complex, -rows)]

fwrite(leocc_cov, file.path(out_dir, "coverage_leocc.csv"))
fwrite(broad_cov, file.path(out_dir, "coverage_broad.csv"))


# ---- Criteria ----------------------------------------------------------------

lc <- rows[file == "LEOCC_Metabolon" & is_complex & !group_name]
bc <- rows[file == "Broad_2022Aug" & is_complex]

crit <- data.table(
  criterion = c("LEOCC complex, non-group rows with one Species ancestor",
                "Broad complex rows parsed to Species level or finer"),
  passed_rows = c(lc[n_species == 1, .N], bc[finer_than_class(name_Level), .N]),
  total_rows  = c(nrow(lc), nrow(bc))
)
crit[, percent := mapply(pct, passed_rows, total_rows)]
crit[, meets_80 := percent >= 80]
fwrite(crit, file.path(out_dir, "criteria.csv"))

# Extra numbers that are not criteria but matter for the design.
extra <- data.table(
  measure = c("LEOCC complex group-name rows (not in the criteria)",
              "LEOCC complex group-name rows whose key maps to one Species",
              "LEOCC complex rows with more than one Species ancestor",
              "LEOCC complex rows whose Species abbreviation rgoslin parsed",
              "LEOCC complex rows where hierarchy and hit abbreviation agree on Species",
              "LEOCC complex rows where they disagree",
              "LEOCC complex rows with no Species that gain exactly one via the connectivity fallback",
              "LEOCC complex rows with one Species after the fallback (original plus fallback)"),
  rows = c(rows[file == "LEOCC_Metabolon" & is_complex & group_name, .N],
           rows[file == "LEOCC_Metabolon" & is_complex & group_name & n_species == 1, .N],
           lc[n_species > 1, .N],
           lc[!is.na(species_from_hierarchy), .N],
           lc[!is.na(species_from_hierarchy) & !is.na(species_from_hit) &
                species_from_hierarchy == species_from_hit, .N],
           lc[!is.na(species_from_hierarchy) & !is.na(species_from_hit) &
                species_from_hierarchy != species_from_hit, .N],
           lc[n_species == 0 & n_species_fallback == 1, .N],
           lc[n_species == 1 | (n_species == 0 & n_species_fallback == 1), .N])
)
fwrite(extra, file.path(out_dir, "extra_measures.csv"))


# ---- Comparison with the manual review ---------------------------------------

review_path <- "notes/review.csv"
if (file.exists(review_path)) {
  review <- fread(review_path, colClasses = "character", na.strings = "")
  labeled <- review[!is.na(category) & nzchar(category)]
  if (nrow(labeled) > 0) {
    long <- labeled[, .(row_id = trimws(unlist(strsplit(row_ids, "|", fixed = TRUE)))),
                    by = .(metlinkr_name, category)]
    check <- merge(long, rows, by = "row_id")
    check <- check[, .(metlinkr_name, category, row_id, name, tier, n_hits,
                       species_from_hierarchy, species_from_fallback, species_from_hit,
                       molecular_from_hit, name_Species.Name,
                       name_Molecular.Species.Name)]
    setorder(check, metlinkr_name, row_id)
    fwrite(check, file.path(out_dir, "review_check.csv"))
    message("Review check: ", uniqueN(check$metlinkr_name),
            " labeled groups contain lipid rows. See review_check.csv")
  }
}


# ---- Print a short summary ---------------------------------------------------

message("\nLipid rows in scope: LEOCC ", rows[file == "LEOCC_Metabolon", .N],
        ", Broad ", rows[file == "Broad_2022Aug", .N])
message("\nCriteria:")
print(crit)
message("\nExtra measures:")
print(extra)
message("\nFiles written to ", out_dir)
