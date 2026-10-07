# Inventory for v0.5: can each row be checked against the RefMet entry
# MetLinkR assigned it to?
#
# MetLinkR's "Standardized name" is a RefMet name, so every group stands
# for one RefMet entry. This script counts, for each row, what evidence
# could test that assignment: the row's own InChIKey against the entry's
# key, or the row's lipid profile against the entry's name parsed with
# rgoslin. It only counts. No expert label is read, and no rule or
# threshold is set here.
#
# Needs: data-raw/refmet/refmet_all.txt (downloaded 2026-10-07), the v0.4.2
# audits of the development data and of COMETS.
# Output: results/refmet_inventory/
# Run from the project root: Rscript scripts/70_refmet_inventory.R

pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(data.table))

refmet_path  <- "data-raw/refmet/refmet_all.txt"
expected_md5 <- "63079d6af6aff5fec9ddc9d3f952afa2"
if (unname(tools::md5sum(refmet_path)) != expected_md5) {
  stop("md5 of ", refmet_path, " differs from the file this inventory was planned on.")
}

out_dir <- "results/refmet_inventory"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- Read the RefMet table ------------------------------------------------
# One "field<TAB>value" line per field, records separated by blank lines,
# each record starting with "name".

x  <- readLines(refmet_path, warn = FALSE, encoding = "UTF-8")
x  <- x[nzchar(x)]
tab <- regexpr("\t", x, fixed = TRUE)
field <- ifelse(tab > 0, substr(x, 1, tab - 1), x)
value <- ifelse(tab > 0, substr(x, tab + 1, nchar(x)), NA_character_)
rec <- cumsum(field == "name")
refmet <- dcast(data.table(rec, field, value), rec ~ field, value.var = "value")
refmet[!is_valid_inchikey(inchi_key), inchi_key := NA_character_]
message("RefMet entries: ", nrow(refmet), ", with a valid InChIKey: ", refmet[!is.na(inchi_key), .N])

dup <- refmet[, .N, by = name][N > 1]
message("Names used by more than one entry: ", nrow(dup))


# ---- Per row: the row's own evidence and its group's RefMet entry ---------

inventory <- function(audit, label) {
  r <- as.data.table(audit$rows)
  need <- c("row_id", "file", "original_name", "group_id", "ik_valid",
            "ik_block1", "inchikey_source", "lipid_species", "name_group",
            "name_alternatives")
  miss <- setdiff(need, names(r))
  if (length(miss) > 0) stop(label, ": audit rows lack ", paste(miss, collapse = ", "))

  r <- r[!is.na(group_id) & group_id != "-"]
  e <- refmet[match(r$group_id, refmet$name)]
  if (anyNA(e$name)) stop(label, ": ", sum(is.na(e$name)), " group names are not in the RefMet table.")
  r[, rm_key := e$inchi_key]
  r[, rm_block1 := substr(rm_key, 1, 14)]

  # The group's RefMet name parsed as a lipid, one call per distinct name.
  u  <- unique(r$group_id)
  lp <- parse_lipid_names(u)
  r[, rm_species := lp$species[match(group_id, u)]]
  r[, rm_level := lp$level[match(group_id, u)]]
  r[!lp$parsed[match(group_id, u)], c("rm_species", "rm_level") := NA_character_]

  r[, row_side := fifelse(name_group | name_alternatives, "group_or_alternative_name",
                  fifelse(ik_valid & !is.na(inchikey_source), paste0("key_", inchikey_source),
                  fifelse(!is.na(lipid_species), "lipid_profile", "none")))]
  r[, check := fifelse(startsWith(row_side, "key_") & !is.na(rm_block1),
                 fifelse(ik_block1 == rm_block1, "key_agrees", "key_differs"),
               fifelse(startsWith(row_side, "key_"), "key_but_entry_has_none",
               fifelse(row_side == "lipid_profile" & !is.na(rm_species),
                 fifelse(lipid_species == rm_species, "species_agrees", "species_differs"),
               fifelse(row_side == "lipid_profile", "lipid_but_entry_not_parsed", "nothing_to_compare"))))]
  r[, dataset := label]
  r
}

dev    <- inventory(readRDS("results/audit_v0.1_sealed.rds"), "development")
comets <- inventory(readRDS("results/comets_v0.4.2/audit.rds"), "comets")
all <- rbind(dev, comets, fill = TRUE)

entry_side <- unique(all[, .(dataset, group_id, has_key = !is.na(rm_key), parsed = !is.na(rm_species))])
message("\nGroup entries: does the RefMet entry carry a key, does its name parse as a lipid?")
print(entry_side[, .N, by = .(dataset, has_key, parsed)][order(dataset, -N)])

message("\nRows by the evidence they carry:")
print(dcast(all[, .N, by = .(dataset, row_side)], row_side ~ dataset, value.var = "N", fill = 0))

message("\nRows by what can be compared with their RefMet entry:")
print(dcast(all[, .N, by = .(dataset, check)], check ~ dataset, value.var = "N", fill = 0))

message("\nCOMETS rows by file and comparison:")
print(dcast(comets[, .N, by = .(file, check)], file ~ check, value.var = "N", fill = 0))

fwrite(all[, .(dataset, row_id, file, original_name, group_id, row_side, check,
               ik_block1, rm_block1, lipid_species, rm_species)],
       file.path(out_dir, "rows.csv"))
fwrite(all[check %in% c("key_differs", "species_differs"),
           .(dataset, file, original_name, group_id, check, inchikey_source,
             ik_block1, rm_block1, lipid_species, rm_species)],
       file.path(out_dir, "disagreements.csv"))
message("\nFiles written to ", out_dir)
