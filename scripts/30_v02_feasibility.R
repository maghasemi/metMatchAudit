# Feasibility test for v0.2: lipid matches across files.
#
# Three questions, counted on the MetLinkR example data (LEOCC Metabolon
# and Broad 2022):
#   1. How many cross-file pairs of lipid rows share a Species? How many
#      of those did MetLinkR already link, and how many are new?
#   2. How often does one row match several rows in the other file?
#   3. For Metabolon group-name rows such as
#      "sphingomyelin (d18:1/20:2, d18:2/20:1)*", do all listed members
#      add up to the same chain sum?
#
# Uses the v0.1 package functions. Nothing here changes them.
# Output: results/v02_feasibility/
# Run from the project root: Rscript scripts/30_v02_feasibility.R

pkgload::load_all(".", quiet = TRUE)

out_dir <- "results/v02_feasibility"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

base <- "data-raw/MetLinkR/inst/extdata/Example_output"
files <- data.frame(
  path         = file.path(base, c("LEOCC_Metabolon_Annotations_metLinkR.csv",
                                   "Broad_2022Aug_annotations_metLinkR.csv")),
  name_col     = c("CHEMICAL_NAME", "Metabolite"),
  inchikey_col = c("INCHIKEY", "inChIKey"),
  label        = c("LEOCC_Metabolon", "Broad_2022Aug"),
  stringsAsFactors = FALSE
)

db   <- load_swisslipids("data-raw/swisslipids/lipids.tsv",
                         expected_md5 = "b6ddd5bc3003967dce95e905e5748953")
rows <- read_metlinkr(files, encoding = "Latin-1")
ev   <- lipid_evidence(row_evidence(rows), db)


# ---- 1. Cross-file pairs with the same Species -------------------------

lip <- ev[!is.na(ev$lipid_species), ]
a <- lip[lip$file == "LEOCC_Metabolon", c("row_id", "lipid_species")]
b <- lip[lip$file == "Broad_2022Aug",   c("row_id", "lipid_species")]
cand <- merge(a, b, by = "lipid_species", suffixes = c("_a", "_b"))

pairs <- data.frame(group_id = cand$lipid_species,
                    row_a = cand$row_id_a, row_b = cand$row_id_b,
                    stringsAsFactors = FALSE)
pairs <- classify_pairs(ev, pairs)

grp <- stats::setNames(ev$group_id, ev$row_id)
pairs$linked_by_metlinkr <- !is.na(grp[pairs$row_a]) & !is.na(grp[pairs$row_b]) &
  grp[pairs$row_a] == grp[pairs$row_b]

rel_table <- as.data.frame(table(relationship = pairs$relationship,
                                 rule = pairs$rule_id,
                                 linked_by_metlinkr = pairs$linked_by_metlinkr))
rel_table <- rel_table[rel_table$Freq > 0, ]


# ---- 2. One-to-one or one-to-many -----------------------------------------

proposed <- pairs[pairs$relationship %in% c("exact", "broader_narrower"), ]
n_a <- table(proposed$row_a)
n_b <- table(proposed$row_b)
proposed$partners_of_a <- as.integer(n_b[proposed$row_b])
proposed$partners_of_b <- as.integer(n_a[proposed$row_a])
proposed$shape <- ifelse(proposed$partners_of_a == 1 & proposed$partners_of_b == 1,
                         "one_to_one", "one_to_many")

nm <- stats::setNames(ev$original_name, ev$row_id)
proposed$name_a <- nm[proposed$row_a]
proposed$name_b <- nm[proposed$row_b]

shape_table <- as.data.frame(table(shape = proposed$shape,
                                   linked_by_metlinkr = proposed$linked_by_metlinkr))


# ---- 3. Group-name lipid rows: do the members share one chain sum? -------
# Only names with a comma-separated list of chain notations inside
# parentheses are looked at, e.g. "(d18:1/20:2, d18:2/20:1, d16:1/22:2)".
# Each member's chains are added up: carbons, double bonds.

chain_sum <- function(member) {
  m <- regmatches(member, gregexpr("[dtmOP]?-?(\\d+):(\\d+)", member))[[1]]
  if (length(m) == 0) return(NA_character_)
  nums <- regmatches(m, regexpr("\\d+:\\d+", m))
  c_ <- sum(as.integer(sub(":.*", "", nums)))
  db_ <- sum(as.integer(sub(".*:", "", nums)))
  paste0(c_, ":", db_)
}

grp_rows <- ev[ev$file == "LEOCC_Metabolon" & ev$name_group, ]
lists <- regmatches(grp_rows$original_name,
                    regexpr("\\([^()]*:[^()]*,[^()]*\\)", grp_rows$original_name))
has_list <- grepl("\\([^()]*:[^()]*,[^()]*\\)", grp_rows$original_name)
grp_rows <- grp_rows[has_list, ]

group_check <- do.call(rbind, lapply(seq_len(nrow(grp_rows)), function(i) {
  inside  <- sub("^\\((.*)\\)$", "\\1", lists[i])
  members <- trimws(strsplit(inside, ",")[[1]])
  sums    <- vapply(members, chain_sum, character(1))
  sl      <- swisslipids_lookup(grp_rows$inchikey[i], db)
  data.frame(row_id = grp_rows$row_id[i],
             name = grp_rows$original_name[i],
             n_members = length(members),
             member_sums = paste(sums, collapse = " | "),
             all_same_sum = length(unique(sums)) == 1 && !anyNA(sums),
             representative_species = sl$sl_species,
             stringsAsFactors = FALSE)
}))


# ---- Write and print ----------------------------------------------------------

utils::write.csv(pairs, file.path(out_dir, "candidate_pairs.csv"), row.names = FALSE)
utils::write.csv(proposed, file.path(out_dir, "proposed_matches.csv"), row.names = FALSE)
utils::write.csv(group_check, file.path(out_dir, "group_rows.csv"), row.names = FALSE)

message("Lipid rows with a profile: LEOCC ", nrow(a), ", Broad ", nrow(b))
message("Cross-file pairs with the same Species: ", nrow(pairs))
message("\n1. Relationship and rule, by whether MetLinkR linked the pair:")
print(rel_table, row.names = FALSE)
message("\n2. Shape of proposed matches (exact or broader_narrower):")
print(shape_table, row.names = FALSE)
message("\n3. Metabolon group-name lipid rows with a chain list: ", nrow(group_check))
message("   all members share one chain sum: ", sum(group_check$all_same_sum))
print(group_check[, c("name", "member_sums", "all_same_sum", "representative_species")],
      row.names = FALSE)
message("\nFiles written to ", out_dir)
