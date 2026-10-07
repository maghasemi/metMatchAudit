# Extracts the official Biocrates p180 old-to-new nomenclature table from
# the vendor's PDF, so the tool never relies on a hand-typed copy.
#
# Source: "biocrates-p180 old vs. new nomenclature (v1-2023)", last
# modified 2023-02-24,
# https://biocrates.com/wp-content/uploads/2023/12/p180-nomenclature.pdf
# The new names follow Liebisch et al. 2020 (LIPID MAPS shorthand).
#
# Output: data-raw/biocrates/p180_nomenclature.csv with one row per
# analyte: section, old (abbreviation as printed), alternative (the
# second abbreviation Biocrates gives in parentheses, if any), new (the
# new name, or NA where the table shows no change).
# The script stops if any section does not have the number of analytes
# the PDF itself states.
#
# Run from the project root: Rscript scripts/10_v0.6/01_biocrates_table.R

pdf_path     <- "data-raw/biocrates/p180-nomenclature-v1-2023.pdf"
expected_md5 <- "4a1f84e43f5ad2f9dfe5aaa093fc0195"
out_path     <- "data-raw/biocrates/p180_nomenclature.csv"

if (unname(tools::md5sum(pdf_path)) != expected_md5) {
  stop("md5 of ", pdf_path, " differs from the file this script was written for.")
}

lines <- unlist(strsplit(pdftools::pdf_text(pdf_path), "\n", fixed = TRUE))

# Section headers as printed. "(continued)" pages belong to the same section.
sections <- c("Amino acids", "Biogenic amines", "Monosaccharides", "Acylcarnitines",
              "Glycerophospholipids", "Sphingolipids")
expected <- c("Amino acids" = 21, "Biogenic amines" = 21, "Monosaccharides" = 1,
              "Acylcarnitines" = 40, "Glycerophospholipids" = 90, "Sphingolipids" = 15)

current <- NA_character_
rows <- list()

for (ln in lines) {
  h <- sub("\\s*\\(continued\\)\\s*$", "", trimws(ln))
  if (h %in% sections) { current <- h; next }
  if (grepl("^biocrates-p180|^For research|Old nomenclature", trimws(ln))) next
  if (is.na(current)) next

  if (current %in% c("Glycerophospholipids", "Sphingolipids")) {
    # e.g. "   PC aa C30:2 2        PC 30:2 2" (the trailing 2 is a footnote)
    m <- regmatches(ln, regexec(
      "^\\s+((?:lysoPC a|PC aa|PC ae|SM \\(OH\\)|SM) C\\d+:\\d+)(?:\\s+\\d)?\\s{2,}((?:LPC|PC O-|PC|SM) ?\\d+:\\d+)(?:\\s+\\d)?\\s*$",
      ln, perl = TRUE))[[1]]
    if (length(m) == 3) rows[[length(rows) + 1]] <- data.frame(
      section = current, old = m[2], alternative = NA_character_,
      new = gsub("\\s+", " ", m[3]), stringsAsFactors = FALSE)
    next
  }

  if (current == "Acylcarnitines") {
    # Main line: "   C5-DC     Glutarylcarnitine", or "   C14:2-OH" followed by dashes.
    # Second line of a pair: "   (C6-OH)   (Hydroxyhexanoylcarnitine)".
    alt <- regmatches(ln, regexec("^\\s+\\((C[0-9][^)]*)\\)", ln))[[1]]
    if (length(alt) == 2 && length(rows) > 0) {
      rows[[length(rows)]]$alternative <- alt[2]
      next
    }
    m <- regmatches(ln, regexec("^\\s+(C\\d+(?::\\d+)?(?:-[A-Z]+)*)\\s", ln))[[1]]
    if (length(m) == 2) rows[[length(rows) + 1]] <- data.frame(
      section = current, old = m[2], alternative = NA_character_, new = NA_character_,
      stringsAsFactors = FALSE)
    next
  }

  # Amino acids, biogenic amines, monosaccharides: abbreviation then name.
  m <- regmatches(ln, regexec("^\\s+(\\S+(?: \\S+)?)\\s{2,}(\\S.*?)\\s{2,}\u2013", ln, perl = TRUE))[[1]]
  if (length(m) == 3) rows[[length(rows) + 1]] <- data.frame(
    section = current, old = sub("\\s+\\d$", "", m[2]), alternative = NA_character_,
    new = NA_character_, stringsAsFactors = FALSE)
}

tab <- do.call(rbind, rows)
counts <- table(factor(tab$section, levels = sections))
print(counts)
bad <- names(expected)[counts[names(expected)] != expected]
if (length(bad) > 0) stop("Analyte counts differ from the PDF for: ", paste(bad, collapse = ", "))
if (anyDuplicated(tab$old)) stop("An old abbreviation appears twice.")

write.csv(tab, out_path, row.names = FALSE)
message("Wrote ", nrow(tab), " analytes to ", out_path)
