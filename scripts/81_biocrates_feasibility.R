# Feasibility for v0.6: can Biocrates p180 rows get a lipid profile from
# the vendor's own nomenclature table, and would that let lipid rules
# reach pairs that now stop at R01?
#
# Counts only. No expert label is read and no rule is changed. Criteria
# are in notes/decisions.md (2026-10-07, v0.6 plan).
#
# Needs: data-raw/biocrates/p180_nomenclature.csv from scripts/80, the
# v0.4.2 COMETS audit.
# Run from the project root: Rscript scripts/81_biocrates_feasibility.R

pkgload::load_all(".", quiet = TRUE)

tab <- read.csv("data-raw/biocrates/p180_nomenclature.csv", colClasses = "character",
                na.strings = c("", "NA"))

# Mapped name per analyte, fixed in the plan:
# - lipids: the new name from the table;
# - acylcarnitines with no alternative given, written Cx, Cx:y or with -OH:
#   CAR x:y, plus ;OH. C0 is carnitine itself and gets nothing. Forms with
#   -DC or -M, and analytes listed with an alternative, get nothing.
ac <- tab$section == "Acylcarnitines"
m  <- regmatches(tab$old, regexec("^C(\\d+)(?::(\\d+))?(-OH)?$", tab$old))
car <- vapply(m, function(v) {
  if (length(v) == 0 || as.integer(v[2]) < 2) return(NA_character_)
  paste0("CAR ", v[2], ":", if (nzchar(v[3])) v[3] else "0", if (nzchar(v[4])) ";OH" else "")
}, character(1))
tab$mapped <- ifelse(tab$section %in% c("Glycerophospholipids", "Sphingolipids"), tab$new,
              ifelse(ac & is.na(tab$alternative), car, NA_character_))

lp <- parse_lipid_names(tab$mapped)
tab$species <- ifelse(lp$parsed, lp$species, NA_character_)
tab$level   <- ifelse(lp$parsed, lp$level, NA_character_)
mapped <- !is.na(tab$mapped)
message("Analytes with a mapped name: ", sum(mapped), " of ", nrow(tab))
message("Mapped names that parse: ", sum(!is.na(tab$species[mapped])), " of ", sum(mapped))
print(table(tab$section[mapped], tab$level[mapped], useNA = "ifany"))
if (any(mapped & is.na(tab$species))) print(tab[mapped & is.na(tab$species), c("old", "mapped")])

# The table key as COMETS writes it: abbreviation, then the alternative in
# parentheses, spaces removed ("C5-DC(C6-OH)").
nospace <- function(x) gsub("\\s+", "", x)
tab$key <- nospace(ifelse(is.na(tab$alternative), tab$old, paste0(tab$old, "(", tab$alternative, ")")))

audit <- readRDS("results/comets_v0.4.2/audit.rds")
rows  <- audit$rows
bio   <- rows$file == "inputs_1539353822"
abbr  <- regmatches(rows$original_name,
                    regexec("\\(\\s*((?:[^()]|\\([^()]*\\))*?)\\s*\\)\\s*$", rows$original_name, perl = TRUE))
rows$bio_key <- NA_character_
rows$bio_key[bio] <- vapply(abbr[bio], function(v) if (length(v) == 2) nospace(v[2]) else NA_character_, character(1))
i <- match(rows$bio_key, tab$key)
rows$bio_section <- tab$section[i]
rows$bio_species <- tab$species[i]

lip_ac <- bio & grepl("carnitine|phosphatidylcholine|sphingomyel", rows$original_name, ignore.case = TRUE)
message("\nCOMETS Biocrates rows: ", sum(bio), "; acylcarnitine and lipid rows: ", sum(lip_ac))
message("Of those, found in the table: ", sum(lip_ac & !is.na(i)),
        sprintf(" (%.1f%%)", 100 * mean(!is.na(i[lip_ac]))))
message("Of those, with a parsed profile: ", sum(lip_ac & !is.na(rows$bio_species)))
message("Not found in the table:")
print(rows$original_name[lip_ac & is.na(i)])

# Cross-file pairs where the Biocrates row would get a profile, and whether
# the other row has a lipid profile today.
sp <- stats::setNames(rows$bio_species, rows$row_id)
ls <- stats::setNames(rows$lipid_species, rows$row_id)
fl <- stats::setNames(rows$file, rows$row_id)
p  <- audit$pairs
p  <- p[fl[p$row_a] != fl[p$row_b] & (!is.na(sp[p$row_a]) | !is.na(sp[p$row_b])), ]
other <- ifelse(!is.na(sp[p$row_a]), p$row_b, p$row_a)
message("\nCross-file pairs with a mapped Biocrates row: ", nrow(p))
message("  other row has a lipid profile: ", sum(!is.na(ls[other])))
message("  other row has no lipid profile: ", sum(is.na(ls[other])))
message("Rule those pairs get today:")
print(table(p$rule_id))

dir.create("results/biocrates_feasibility", showWarnings = FALSE, recursive = TRUE)
write.csv(tab, "results/biocrates_feasibility/table_mapped.csv", row.names = FALSE)
