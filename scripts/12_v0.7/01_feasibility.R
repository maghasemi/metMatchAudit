# v0.7 feasibility (notes/decisions.md, 2026-10-07, v0.7 plan, step 2).
# Changes nothing in the package. Counts, on the development data, COMETS
# and the two files of the second blind test:
#   1. names carrying an isomer marker,
#   2. how many of B's lipid rows get a profile from their name once the
#      name is cleaned (the 80% threshold for change 2),
#   3. names with an optical rotation prefix.
#
# Reads the v0.6.0 audits already on disk. Lists that name rows of the
# private file A go to data-raw/own/blind2/v07_feasibility/; everything
# else to results/v07_feasibility/.
#
# Run from the project root: Rscript scripts/12_v0.7/01_feasibility.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.6.0")
options(width = 200)

pub  <- "results/v07_feasibility"
priv <- "data-raw/own/blind2/v07_feasibility"
dir.create(pub,  showWarnings = FALSE, recursive = TRUE)
dir.create(priv, showWarnings = FALSE, recursive = TRUE)

dev  <- readRDS("results/v06_check/audit_development.rds")$rows
com  <- readRDS("results/v06_check/audit_comets.rds")$rows
bl   <- readRDS("data-raw/own/blind2/results/audit.rds")$rows
sets <- list(development = dev, comets = com,
             blind2_A = bl[bl$file == "A", ], blind2_B = bl[bl$file == "B", ])


# ---- Prototypes of the three changes ---------------------------------------

# 1. Isomer marker: the word "isomer" or "isomers". "stereoisomer" does not
#    match, since the word must start at a boundary.
has_isomer_marker <- function(x) grepl("\\bisomers?\\b", x, ignore.case = TRUE, perl = TRUE)

# 3. Optical rotation prefix.
has_rotation_prefix <- function(x) grepl("^\\s*\\((\\+|-|±|\\+/-)\\)-", x, perl = TRUE)

# 2. MS-DIAL name cleaning. Returns the annotation to parse, or NA for
#    names that list alternatives.
clean_msdial <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("[;\\s]+$", "", x, perl = TRUE)
  alt <- grepl("\\sor\\s", x, perl = TRUE)
  parts <- strsplit(x, "\\s*;\\s*", perl = TRUE)
  last <- vapply(parts, function(p) { p <- p[nzchar(p)]; if (length(p) == 0) NA_character_ else p[length(p)] },
                 character(1))
  last <- sub("^ACar\\s", "CAR ", last, perl = TRUE)
  ifelse(alt, NA_character_, last)
}


# ---- 1 and 3: counts in all four sets --------------------------------------

marker <- do.call(rbind, lapply(names(sets), function(s) {
  r <- sets[[s]]
  data.frame(set = s, rows = nrow(r),
             isomer_marker = sum(has_isomer_marker(r$original_name)),
             isomer_marker_with_key = sum(has_isomer_marker(r$original_name) & r$ik_valid),
             rotation_prefix = sum(has_rotation_prefix(r$original_name)),
             stringsAsFactors = FALSE)
}))
message("==== Isomer markers and rotation prefixes")
print(marker, row.names = FALSE, right = FALSE)
write.csv(marker, file.path(pub, "marker_counts.csv"), row.names = FALSE)

listing <- function(test) do.call(rbind, lapply(names(sets), function(s) {
  r <- sets[[s]]; k <- test(r$original_name)
  if (!any(k)) return(NULL)
  data.frame(set = s, row_id = r$row_id[k], name = r$original_name[k],
             inchikey = r$inchikey[k], group_id = r$group_id[k], stringsAsFactors = FALSE)
}))
for (what in c("isomer", "rotation")) {
  l <- listing(if (what == "isomer") has_isomer_marker else has_rotation_prefix)
  if (is.null(l)) next
  message("\n", what, " names outside file A:")
  print(l[l$set != "blind2_A", c("set", "row_id", "name")], row.names = FALSE, right = FALSE)
  message(what, " names in file A: ", sum(l$set == "blind2_A"), " (listed in ", priv, ")")
  write.csv(l[l$set != "blind2_A", ], file.path(pub, paste0(what, "_names.csv")), row.names = FALSE, na = "")
  write.csv(l[l$set == "blind2_A", ], file.path(priv, paste0(what, "_names_A.csv")), row.names = FALSE, na = "")
}


# ---- 2: MS-DIAL names in B -------------------------------------------------

lip_pat <- "^(L?P[ACEGIS]|SM|Cer|HexCer|Hex2Cer|GlcCer|TAG|TG|DAG|DG|MAG|MG|CE|ACar|CAR|FA)[ -]"
b <- sets$blind2_B
b <- b[grepl(lip_pat, trimws(b$original_name)), ]
cleaned <- clean_msdial(b$original_name)
p <- parse_lipid_names(cleaned)

key_source <- b$lipid_source %in% c("swisslipids", "swisslipids_fallback", "lmsd")
agree <- p$parsed & key_source & p$species == b$lipid_species
clash <- p$parsed & key_source & p$species != b$lipid_species

res <- data.frame(
  measure = c("B lipid rows", "rows flagged as group or alternative names now",
              "names with alternatives (' or '), no profile by design",
              "cleaned names that parse", "  share of rows that could parse",
              "parse and key profile agree on Species", "parse and key profile disagree (would lose the profile)",
              "parse where the row has no profile now"),
  value = c(nrow(b), sum(b$name_group | b$name_alternatives),
            sum(is.na(cleaned)),
            sum(p$parsed), round(100 * sum(p$parsed) / sum(!is.na(cleaned)), 1),
            sum(agree), sum(clash), sum(p$parsed & is.na(b$lipid_species))),
  stringsAsFactors = FALSE)
message("\n==== MS-DIAL names in B (threshold: 80% of rows that could parse)")
print(res, row.names = FALSE, right = FALSE)
message("\nLevels from the cleaned names:"); print(table(p$level, useNA = "ifany"))
message("Grammar of names that did not parse:"); print(table(p$grammar[!p$parsed & !is.na(cleaned)], useNA = "ifany"))

fail <- data.frame(row_id = b$row_id, name = b$original_name, cleaned = cleaned,
                   grammar = p$grammar, stringsAsFactors = FALSE)[!p$parsed & !is.na(cleaned), ]
fail$class <- sub("[ -].*$", "", trimws(fail$name))
message("\nNames that did not parse, by class:"); print(sort(table(fail$class), decreasing = TRUE))
message("\nFirst 30 that did not parse:")
print(utils::head(fail[, c("row_id", "name", "cleaned")], 30), row.names = FALSE, right = FALSE)

message("\nDisagreements between the cleaned name and the key profile:")
dis <- data.frame(row_id = b$row_id, name = b$original_name, name_species = p$species,
                  key_species = b$lipid_species, key_source = b$lipid_source, stringsAsFactors = FALSE)[clash, ]
print(utils::head(dis, 30), row.names = FALSE, right = FALSE)

write.csv(res, file.path(pub, "msdial_counts.csv"), row.names = FALSE)
write.csv(fail, file.path(pub, "msdial_not_parsed.csv"), row.names = FALSE, na = "")
write.csv(dis,  file.path(pub, "msdial_species_disagree.csv"), row.names = FALSE, na = "")
message("\nCounts and B lists in ", pub, "; A lists in ", priv, ".")
