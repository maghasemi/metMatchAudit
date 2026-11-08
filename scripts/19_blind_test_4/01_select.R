# Fourth blind test, step 2: chooses and copies the three files under the
# rule in notes/decisions.md (2026-11-08, fourth blind test, plan).
#
#   1. Candidates: studies the Metabolomics Workbench lists for MS, human,
#      blood.
#   2. An analysis qualifies with at least 100 named metabolites.
#   3. It is left out when more than half of its names equal the Workbench's
#      RefMet name after normalize_name().
#   4. Three analyses from three institutes and three studies: one whose
#      table gives a PubChem CID, InChIKey or HMDB ID for at least half of
#      its rows, and two with no such column.
#   5. Only the lab's name and its identifier columns are kept.
#
# Prints no metabolite name. Every answer from the Workbench is cached in
# data-raw/blind4/cache/, so a rerun gives the same choice. The copies go
# to data-raw/blind4/ (git ignores data-raw/); the log of the choice, which
# holds no names, goes to results/blind4_selection/.
#
# Needs internet. Expect several minutes: a few hundred studies are asked.
#
# Run from the project root:
#   Rscript scripts/19_blind_test_4/01_select.R 2>&1 | tee results/blind4_selection_log.txt

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.11.0")
options(width = 200)

raw <- "data-raw/blind4"
cache_dir <- file.path(raw, "cache")
out <- "results/blind4_selection"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (length(list.files(raw, pattern = "^F[123]\\.csv$")) > 0) {
  stop("The files were chosen already. Delete data-raw/blind4/F*.csv only on purpose.")
}

base <- "https://www.metabolomicsworkbench.org/rest/"
n_failed <- 0

# One request, cached on disk by its path. A failed request is not cached
# and gives NULL.
ask <- function(path) {
  f <- file.path(cache_dir, paste0(gsub("[^A-Za-z0-9]+", "_", path), ".rds"))
  if (file.exists(f)) return(readRDS(f))
  resp <- tryCatch({
    req <- httr2::request(paste0(base, path))
    req <- httr2::req_user_agent(req, "metMatchAudit")
    req <- httr2::req_timeout(req, 120)
    req <- httr2::req_retry(req, max_tries = 3)
    httr2::resp_body_string(httr2::req_perform(req))
  }, error = function(e) NULL)
  Sys.sleep(0.2)
  if (is.null(resp)) {
    n_failed <<- n_failed + 1
    return(NULL)
  }
  saveRDS(resp, f)
  resp
}
ask_json <- function(path) {
  txt <- ask(path)
  if (is.null(txt) || !nzchar(trimws(txt)) || trimws(txt) %in% c("[]", "{}")) return(NULL)
  tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
}

# The Workbench returns one record as an object and several as an object of
# numbered records. Both become a list of records.
records <- function(x) {
  if (is.null(x)) return(list())
  if (!is.null(x$study_id) || !is.null(x$study)) list(x) else unname(x)
}
field <- function(r, name) if (is.null(r[[name]])) NA_character_ else as.character(r[[name]])


# ---- 1. Candidate studies ------------------------------------------------------

studies <- unique(vapply(records(ask_json("metstat/MS;;;Human;Blood;;;")), field, "", "study"))
studies <- sort(studies[!is.na(studies)])
message("1. Studies listed for MS, human, blood: ", length(studies))


# ---- 2 and 3. Analyses with enough names, not already RefMet names ------------

# Share of names that are already the RefMet name. A missing RefMet name
# counts as different.
same_as_refmet <- function(name, refmet) {
  a <- normalize_name(name)
  b <- normalize_name(refmet)
  mean(!is.na(a) & !is.na(b) & a == b)
}

rows <- list()
for (i in seq_along(studies)) {
  s <- studies[i]
  sm <- records(ask_json(paste0("study/study_id/", s, "/summary")))
  info <- if (length(sm) > 0) sm[[1]] else list()
  met <- records(ask_json(paste0("study/study_id/", s, "/metabolites")))
  if (length(met) > 0) {
    m <- data.frame(analysis_id = vapply(met, field, "", "analysis_id"),
                    name = vapply(met, field, "", "metabolite_name"),
                    refmet = vapply(met, field, "", "refmet_name"), stringsAsFactors = FALSE)
    m <- m[!is.na(m$name) & nzchar(trimws(m$name)), ]
    for (a in unique(m$analysis_id)) {
      x <- m[m$analysis_id == a, ]
      rows[[length(rows) + 1]] <- data.frame(
        study = s, analysis_id = a, institute = field(info, "institute"),
        licence = field(info, "license"), release_date = field(info, "release_date"),
        title = field(info, "study_title"), n_named = nrow(x),
        refmet_share = same_as_refmet(x$name, x$refmet),
        stringsAsFactors = FALSE)
    }
  }
  if (i %% 50 == 0) message("  ", i, " of ", length(studies), " studies")
}
cand <- do.call(rbind, rows)
message("Analyses with at least one named metabolite: ", nrow(cand))
cand$enough <- cand$n_named >= 100
cand$refmet_named <- cand$refmet_share > 0.5
message("2. With at least 100 named metabolites: ", sum(cand$enough))
message("3. Of those, left out because more than half the names are RefMet names: ",
        sum(cand$enough & cand$refmet_named))
q <- cand[cand$enough & !cand$refmet_named, ]


# ---- 4. Identifier columns -----------------------------------------------------

# The METABOLITES block of the mwTab file holds the lab's names and any
# identifier columns it chose to give.
metabolite_table <- function(a) {
  txt <- ask(paste0("study/analysis_id/", a, "/mwtab/txt"))
  if (is.null(txt)) return(NULL)
  lines <- strsplit(txt, "\r?\n")[[1]]
  s <- grep("^METABOLITES_START", lines)
  e <- grep("^METABOLITES_END", lines)
  if (length(s) != 1 || length(e) != 1 || e - s < 3) return(NULL)
  body <- strsplit(lines[(s + 1):(e - 1)], "\t", fixed = TRUE)
  head <- trimws(body[[1]])
  width <- length(head)
  cells <- lapply(body[-1], function(v) { v <- trimws(v); length(v) <- width; v })
  d <- as.data.frame(do.call(rbind, cells), stringsAsFactors = FALSE)
  names(d) <- make.unique(head)
  d
}
id_columns <- function(d) {
  h <- tolower(names(d))
  list(pubchem = names(d)[grepl("pubchem|^cid$", h)][1],
       inchikey = names(d)[grepl("inchi.?key", h)][1],
       hmdb = names(d)[grepl("hmdb", h)][1])
}
id_coverage <- function(d) {
  cols <- id_columns(d)
  ok <- rep(FALSE, nrow(d))
  if (!is.na(cols$pubchem)) ok <- ok | !is.na(normalize_pubchem_id(d[[cols$pubchem]]))
  if (!is.na(cols$inchikey)) ok <- ok | is_valid_inchikey(d[[cols$inchikey]])
  if (!is.na(cols$hmdb)) ok <- ok | !is.na(normalize_hmdb_id(d[[cols$hmdb]]))
  c(has_column = any(!is.na(unlist(cols))), coverage = if (nrow(d) > 0) mean(ok) else 0)
}

message("\nReading the metabolite tables of ", nrow(q), " analyses")
q$has_id_column <- NA
q$id_coverage <- NA_real_
for (i in seq_len(nrow(q))) {
  d <- metabolite_table(q$analysis_id[i])
  if (!is.null(d)) {
    cv <- id_coverage(d)
    q$has_id_column[i] <- as.logical(cv[["has_column"]])
    q$id_coverage[i] <- cv[["coverage"]]
  }
  if (i %% 50 == 0) message("  ", i, " of ", nrow(q))
}
q$stratum <- ifelse(is.na(q$has_id_column), "no table",
             ifelse(q$has_id_column & q$id_coverage >= 0.5, "identifiers",
             ifelse(!q$has_id_column, "names only", "some identifiers")))
message("4. Qualifying analyses by stratum:")
print(table(q$stratum))
if (n_failed > 0) stop(n_failed, " requests failed. Run the script again before anything is drawn.")


# ---- Draw ------------------------------------------------------------------------

# A study is drawn first and then one of its analyses, so a study with many
# analyses is not more likely to be chosen.
inst_key <- function(x) normalize_name(gsub("[[:punct:]]", " ", x))
set.seed(20261108)
pick <- function(pool, taken) {
  pool <- pool[!pool$study %in% taken$study & !inst_key(pool$institute) %in% inst_key(taken$institute) &
                 !is.na(pool$institute), ]
  if (nrow(pool) == 0) stop("No analysis left in this stratum under the rule.")
  s <- sample(unique(pool$study), 1)
  a <- pool[pool$study == s, ]
  a[sample(nrow(a), 1), ]
}
chosen <- q[0, ]
chosen <- rbind(chosen, pick(q[q$stratum == "identifiers", ], chosen))
chosen <- rbind(chosen, pick(q[q$stratum == "names only", ], chosen))
chosen <- rbind(chosen, pick(q[q$stratum == "names only", ], chosen))
chosen$label <- c("F1", "F2", "F3")
chosen$rows_copied <- NA_integer_
chosen$columns <- NA_character_
chosen$md5 <- NA_character_
q$chosen <- q$analysis_id %in% chosen$analysis_id


# ---- 5. Copies -------------------------------------------------------------------

for (i in seq_len(nrow(chosen))) {
  d <- metabolite_table(chosen$analysis_id[i])
  cols <- id_columns(d)
  name_col <- names(d)[1]
  copy <- data.frame(BIOCHEMICAL = d[[name_col]], stringsAsFactors = FALSE)
  if (!is.na(cols$pubchem)) copy$PUBCHEM <- d[[cols$pubchem]]
  if (!is.na(cols$hmdb)) copy$HMDB <- d[[cols$hmdb]]
  if (!is.na(cols$inchikey)) copy$INCHIKEY <- d[[cols$inchikey]]
  copy <- copy[!is.na(copy$BIOCHEMICAL) & nzchar(copy$BIOCHEMICAL), , drop = FALSE]
  path <- file.path(raw, paste0(chosen$label[i], ".csv"))
  write.csv(copy, path, row.names = FALSE, na = "")
  chosen$rows_copied[i] <- nrow(copy)
  chosen$columns[i] <- paste(names(copy), collapse = " ")
  chosen$md5[i] <- unname(tools::md5sum(path))
}

message("\n==== chosen")
print(chosen[, c("label", "stratum", "study", "analysis_id", "institute", "licence", "release_date",
                 "n_named", "refmet_share", "id_coverage", "rows_copied", "columns", "md5")],
      row.names = FALSE, right = FALSE)
message("Titles:")
for (i in seq_len(nrow(chosen))) message("  ", chosen$label[i], ": ", chosen$title[i])

write.csv(cand, file.path(out, "analyses.csv"), row.names = FALSE)
write.csv(q, file.path(out, "qualifying.csv"), row.names = FALSE)
write.csv(chosen, file.path(out, "chosen.csv"), row.names = FALSE)
message("\nFailed requests: ", n_failed)
