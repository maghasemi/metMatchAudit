# Structure evidence through RDKit (optional).
#
# RDKit is reached through reticulate. It is used only for the v0.3
# checks: rule R09b and row check C04. Everything else in the package
# works without it.

.rdkit_env <- new.env(parent = emptyenv())

# Returns the RDKit modules, or NULL when RDKit cannot be loaded. The
# result is cached for the session.
.rdkit <- function() {
  if (exists("mods", envir = .rdkit_env)) return(get("mods", envir = .rdkit_env))
  mods <- NULL
  if (requireNamespace("reticulate", quietly = TRUE)) {
    ns <- asNamespace("reticulate")
    if (exists("py_require", envir = ns)) {
      try(get("py_require", envir = ns)("rdkit"), silent = TRUE)
    }
    mods <- tryCatch({
      m <- list(Chem  = reticulate::import("rdkit.Chem", delay_load = FALSE),
                Inchi = reticulate::import("rdkit.Chem.inchi", delay_load = FALSE))
      reticulate::import("rdkit.RDLogger")$DisableLog("rdApp.*")
      m
    }, error = function(e) NULL)
  }
  assign("mods", mods, envir = .rdkit_env)
  mods
}

#' Is RDKit available?
#'
#' @return `TRUE` if RDKit can be loaded through reticulate.
#' @export
rdkit_available <- function() !is.null(.rdkit())

.smiles_to_inchi <- function(smiles, rd) {
  vapply(smiles, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    # Drop the CXSMILES extension after a space; RDKit would read it,
    # but relative stereo is not what the stored keys describe.
    mol <- rd$Chem$MolFromSmiles(sub(" .*$", "", s))
    if (is.null(mol)) return(NA_character_)
    out <- rd$Inchi$MolToInchi(mol)
    if (is.null(out) || !nzchar(out)) NA_character_ else out
  }, character(1), USE.NAMES = FALSE)
}

.inchi_to_key <- function(inchi, rd) {
  vapply(inchi, function(x) {
    if (is.na(x) || !startsWith(x, "InChI=")) return(NA_character_)
    k <- tryCatch(rd$Inchi$InchiToInchiKey(x), error = function(e) NULL)
    if (is.null(k) || !nzchar(k)) NA_character_ else k
  }, character(1), USE.NAMES = FALSE)
}

#' Add structure evidence from SMILES or InChI
#'
#' For each row, takes the InChI from an `inchi` column if there is one,
#' otherwise builds it from a `smiles` column with RDKit. The InChIKey
#' computed from that structure is compared with the row's stored key.
#' The structure is trusted (`st_usable_inchi`) only when the two agree
#' on their first 25 characters, so the structure and the key describe
#' the same molecule.
#'
#' Needs RDKit through reticulate. Without it, the columns are added but
#' left empty, and a message says so.
#'
#' @param evidence Rows after [row_evidence()].
#' @param use_rdkit Set to `FALSE` to skip this step.
#'
#' @return `evidence` with columns `st_inchi`, `st_source` (`"file_inchi"`
#'   or `"rdkit_smiles"`), `st_key`, `st_key_agrees`, `st_block1_differs`
#'   and `st_usable_inchi`.
#'
#' @export
structure_evidence <- function(evidence, use_rdkit = TRUE) {
  new_cols <- c("st_inchi", "st_source", "st_key", "st_key_agrees",
                "st_block1_differs", "st_usable_inchi")
  clash <- intersect(new_cols, names(evidence))
  if (length(clash) > 0) {
    stop("`evidence` already has these columns: ", paste(clash, collapse = ", "), ".",
         call. = FALSE)
  }
  n <- nrow(evidence)
  evidence$st_inchi <- NA_character_
  evidence$st_source <- NA_character_
  evidence$st_key <- NA_character_
  evidence$st_key_agrees <- NA
  evidence$st_block1_differs <- NA
  evidence$st_usable_inchi <- NA_character_
  if (!use_rdkit || n == 0) return(evidence)

  rd <- .rdkit()
  if (is.null(rd)) {
    message("RDKit is not available, so structure checks (R09b, C04) are skipped. ",
            "Install it with reticulate::py_require(\"rdkit\").")
    return(evidence)
  }

  file_inchi <- if ("inchi" %in% names(evidence)) as.character(evidence$inchi) else rep(NA_character_, n)
  file_inchi[!is.na(file_inchi) & !startsWith(file_inchi, "InChI=")] <- NA_character_
  smiles <- if ("smiles" %in% names(evidence)) as.character(evidence$smiles) else rep(NA_character_, n)

  use_file <- !is.na(file_inchi)
  inchi <- file_inchi
  need_smiles <- !use_file & !is.na(smiles)
  inchi[need_smiles] <- .smiles_to_inchi(smiles[need_smiles], rd)

  key <- .inchi_to_key(inchi, rd)
  stored <- as.character(evidence$inchikey)
  comparable <- !is.na(key) & is_valid_inchikey(stored)

  evidence$st_inchi <- inchi
  evidence$st_source <- ifelse(use_file, "file_inchi",
                               ifelse(!is.na(inchi), "rdkit_smiles", NA_character_))
  evidence$st_key <- key
  evidence$st_key_agrees <- ifelse(comparable, substr(key, 1, 25) == substr(stored, 1, 25), NA)
  evidence$st_block1_differs <- ifelse(comparable, substr(key, 1, 14) != substr(stored, 1, 14), NA)
  evidence$st_usable_inchi <- ifelse(comparable & evidence$st_key_agrees, inchi, NA_character_)
  evidence
}
