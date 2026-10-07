# Lipid profiles from a vendor's own nomenclature table. So far only the
# Biocrates AbsoluteIDQ p180 kit (notes/decisions.md, 2026-10-07, v0.6
# plan). The table ships in inst/extdata and was extracted from the
# vendor's PDF by scripts/80_biocrates_table.R.

.known_vendors <- c("biocrates_p180")

#' The Biocrates p180 nomenclature table
#'
#' The old-to-new nomenclature table Biocrates published for the
#' AbsoluteIDQ p180 kit ("biocrates-p180 old vs. new nomenclature
#' (v1-2023)", last modified 2023-02-24). Its new lipid names follow the
#' LIPID MAPS shorthand of Liebisch et al. 2020.
#'
#' Each analyte gets a mapped name:
#'
#' * glycerophospholipids and sphingolipids: the new name from the table,
#'   as printed (so `lysoPC a C28:0` becomes `PC O-28:0` and
#'   `SM (OH) C22:1` becomes `SM 41:1`, as Biocrates decided);
#' * acylcarnitines written `Cx`, `Cx:y` or `Cx:y-OH`, with x of at least
#'   2 and no alternative listed: `CAR x:y`, plus `;OH`;
#' * everything else (carnitine, the -DC and -M forms, analytes listed
#'   with an alternative, amino acids, amines, hexoses): nothing.
#'
#' @return A data frame with `section`, `old`, `alternative`, `new`,
#'   `mapped`, `species` and `key`, the form a row name is matched on:
#'   the abbreviation with its alternative in parentheses, spaces
#'   removed.
#'
#' @export
biocrates_p180_table <- function() {
  path <- system.file("extdata", "biocrates_p180_nomenclature.csv", package = "metMatchAudit")
  if (!nzchar(path)) stop("The Biocrates table is missing from the package.", call. = FALSE)
  tab <- utils::read.csv(path, colClasses = "character", na.strings = c("", "NA"))

  m <- regmatches(tab$old, regexec("^C(\\d+)(?::(\\d+))?(-OH)?$", tab$old))
  car <- vapply(m, function(v) {
    if (length(v) == 0 || as.integer(v[2]) < 2) return(NA_character_)
    paste0("CAR ", v[2], ":", if (nzchar(v[3])) v[3] else "0", if (nzchar(v[4])) ";OH" else "")
  }, character(1))
  lipid <- tab$section %in% c("Glycerophospholipids", "Sphingolipids")
  ac <- tab$section == "Acylcarnitines" & is.na(tab$alternative)
  tab$mapped <- ifelse(lipid, tab$new, ifelse(ac, car, NA_character_))

  lp <- parse_lipid_names(tab$mapped)
  tab$species <- ifelse(lp$parsed, lp$species, NA_character_)
  tab$key <- gsub("\\s+", "", ifelse(is.na(tab$alternative), tab$old,
                                     paste0(tab$old, "(", tab$alternative, ")")))
  tab
}

# The abbreviation a vendor row carries in its last pair of parentheses,
# as in "Tetradecadienylcarnitine  ( C14:2 )" or
# "Glutarylcarnitine (Hydroxyhexanoylcarnitine)  ( C5-DC(C6-OH) )".
# A name without parentheses is taken whole.
.vendor_abbreviation <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexec("\\(\\s*((?:[^()]|\\([^()]*\\))*?)\\s*\\)\\s*$", x, perl = TRUE))
  out <- vapply(seq_along(x), function(i) {
    if (length(m[[i]]) == 2) m[[i]][2] else x[i]
  }, character(1))
  gsub("\\s+", "", out)
}

#' Lipid profiles from a vendor nomenclature table
#'
#' For rows from a vendor file, looks up the row's abbreviation in that
#' vendor's table (see [biocrates_p180_table()]) and sets the lipid
#' profile from the mapped name. The level is always `"species"`:
#' Biocrates measures these analytes by flow injection and cannot tell
#' isomers apart, so the value is a sum whatever rgoslin reports for the
#' name. Only analytes in the table are mapped.
#'
#' A row that already has a profile from another source keeps it when
#' the Species agree. When they differ, the row loses its profile and is
#' flagged in `lipid_conflict`, as with SwissLipids against the name.
#'
#' @param evidence Rows after [lipid_evidence()].
#' @param vendor A character vector, one value per row: `"biocrates_p180"`
#'   or `NA`.
#'
#' @return `evidence` with the profile columns updated and `lipid_source`
#'   set to `"vendor_table"` where the table was used.
#'
#' @export
vendor_lipid_profiles <- function(evidence, vendor) {
  if (length(vendor) != nrow(evidence)) stop("`vendor` needs one value per row.", call. = FALSE)
  bad <- setdiff(stats::na.omit(unique(vendor)), .known_vendors)
  if (length(bad) > 0) {
    stop("Unknown vendor: ", paste(bad, collapse = ", "), ". Known: ",
         paste(.known_vendors, collapse = ", "), ".", call. = FALSE)
  }
  use <- which(!is.na(vendor) & vendor == "biocrates_p180")
  if (length(use) == 0) return(evidence)

  tab <- biocrates_p180_table()
  sp <- tab$species[match(.vendor_abbreviation(evidence$original_name[use]), tab$key)]
  has <- !is.na(sp)
  i <- use[has]; sp <- sp[has]

  old <- evidence$lipid_species[i]
  clash <- !is.na(old) & old != sp
  set <- !clash & is.na(old)
  evidence$lipid_conflict[i[clash]]  <- TRUE
  evidence$lipid_species[i[clash]]   <- NA_character_
  evidence$lipid_molecular[i[clash]] <- NA_character_
  evidence$lipid_level[i[clash]]     <- NA_character_
  evidence$lipid_source[i[clash]]    <- NA_character_

  evidence$lipid_species[i[set]]   <- sp[set]
  evidence$lipid_molecular[i[set]] <- NA_character_
  evidence$lipid_level[i[set]]     <- "species"
  evidence$lipid_source[i[set]]    <- "vendor_table"
  evidence
}
