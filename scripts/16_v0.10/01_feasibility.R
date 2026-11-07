# v0.10 feasibility (notes/decisions.md, 2026-11-07, v0.10 plan, step 2).
# Changes nothing in the package. Asks PubChem for row names and works out
# what a key from the name would do to the R01 pairs.
#
# Only the development data and COMETS are used. Names from the blind test
# files are never sent to PubChem.
#
# Needs internet. About 2,000 names, two requests each, so expect about 20
# minutes. Answers are cached, so an interrupted run can simply be started
# again, and a second run takes seconds.
#
# Run from the project root:
#   Rscript scripts/16_v0.10/01_feasibility.R 2>&1 | tee results/v10_feasibility_log.txt

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.9.0")
options(width = 200)

out <- "results/v10_feasibility"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

dev <- readRDS("results/v09_check/audit_development.rds")
com <- readRDS("results/v09_check/audit_comets_no_fallback.rds")

fetch <- metMatchAudit:::.pubchem_fetch
cache <- metMatchAudit:::.new_cache(file.path(out, "pubchem_name_cache.rds"))


# ---- Which rows --------------------------------------------------------------

clean_name <- function(r) {
  !is.na(r$original_name) & nzchar(trimws(r$original_name)) & !r$name_group & !r$name_alternatives &
    !r$name_internal_standard & !r$name_isomer_marker & !r$name_isomer_label & !r$name_tentative
}
pick <- function(audit, set) {
  r <- audit$rows
  ok <- clean_name(r) & is.na(r$lipid_species)
  data.frame(set = set, row_id = r$row_id[ok], name = r$original_name[ok], ik_valid = r$ik_valid[ok],
             inchikey = r$inchikey[ok], inchikey_source = r$inchikey_source[ok],
             states_stereo = r$name_states_stereo[ok], stringsAsFactors = FALSE)
}
rows <- rbind(pick(dev, "development"), pick(com, "comets"))
rows$role <- ifelse(rows$ik_valid, "check", "fill")
message("Rows with a clean name and no lipid profile:")
print(table(set = rows$set, role = rows$role))
message("Keys of the rows used as a check, by source:")
print(table(set = rows$set[rows$role == "check"], source = rows$inchikey_source[rows$role == "check"], useNA = "ifany"))


# ---- Ask PubChem -------------------------------------------------------------

# A not-found answer is kept. A busy server or a network error is not, so the
# next run asks again.
ask_name <- function(name) {
  key <- paste("name:", name)
  hit <- cache$env[[key]]
  if (!is.null(hit)) return(hit)
  res <- fetch("compound/name/cids/JSON", list(name = name))
  code <- res$Fault$Code
  if (!is.null(res$error) || (!is.null(code) && !identical(code, "PUGREST.NotFound"))) {
    return(list(status = "error", cids = character(0)))
  }
  cids <- as.character(res$IdentifierList$CID)
  ans <- list(status = if (length(cids) > 0) "found" else "not_found", cids = cids)
  cache$env[[key]] <- ans
  saveRDS(cache$env, cache$file)
  ans
}

# One word in a last bracket, at least two characters: "( Pro-PTC )", "(AMP)".
short_name <- function(x) trimws(sub("\\s*\\(\\s*[A-Za-z][A-Za-z0-9-]+\\s*\\)\\s*$", "", x))

names_all <- unique(rows$name)
message("\nAsking PubChem for ", length(names_all), " names")
answers <- vector("list", length(names_all))
names(answers) <- names_all
for (i in seq_along(names_all)) {
  nm <- names_all[i]
  a <- ask_name(nm)
  a$sent <- nm
  stereo <- isTRUE(any(rows$states_stereo[rows$name == nm]))
  shorter <- short_name(nm)
  if (a$status == "not_found" && !stereo && nzchar(shorter) && shorter != nm) {
    a <- ask_name(shorter)
    a$sent <- shorter
  }
  answers[[i]] <- a
  if (i %% 200 == 0) message("  ", i, " of ", length(names_all))
}

# Parent of every CID, then the keys of all parents in batches.
all_cids <- unique(unlist(lapply(answers, `[[`, "cids")))
message("Finding parents of ", length(all_cids), " CIDs")
parents <- vapply(all_cids, metMatchAudit:::.parent_cid, character(1), cache = cache, fetch = fetch)
keys <- metMatchAudit:::.keys_for_cids(unique(parents), cache, fetch)

one_answer <- function(a) {
  if (a$status != "found") return(c(status = a$status, key = NA_character_))
  k <- unique(stats::na.omit(unname(keys[parents[a$cids]])))
  if (length(k) == 1) c(status = "one_key", key = k)
  else c(status = if (length(k) == 0) "no_key" else "several_keys", key = NA_character_)
}
res <- t(vapply(answers, one_answer, c(status = "", key = "")))
i <- match(rows$name, names_all)
rows$status <- res[i, "status"]
rows$name_key <- res[i, "key"]
rows$sent <- vapply(answers, `[[`, "", "sent")[i]
rows$shortened <- rows$sent != rows$name
n_err <- sum(rows$status == "error")
if (n_err > 0) message("Rows with a server or network error: ", n_err, ". Run the script again before reading the results.")


# ---- 1. Rows that already have a key -----------------------------------------

block1 <- function(k) substr(k, 1, 14)
chk <- rows[rows$role == "check", ]
ans <- chk[chk$status == "one_key", ]
ans$same_block1 <- block1(ans$name_key) == block1(ans$inchikey)
ans$same_key <- ans$name_key == ans$inchikey
share_answered <- nrow(ans) / nrow(chk)
share_block1 <- mean(ans$same_block1)

message("\n1. Rows with a key: PubChem gave one key for ", nrow(ans), " of ", nrow(chk),
        " (", round(100 * share_answered, 1), "%). Needed: at least 50%.")
print(table(set = chk$set, status = chk$status))
message("Same connectivity: ", sum(ans$same_block1), " of ", nrow(ans), " (", round(100 * share_block1, 1),
        "%). Needed: at least 98%.")
message("Identical key: ", sum(ans$same_key), " of ", nrow(ans), " (", round(100 * mean(ans$same_key), 1), "%).")
print(table(set = ans$set, same_connectivity = ans$same_block1))
bad <- ans[!ans$same_block1, c("set", "name", "sent", "inchikey", "name_key", "inchikey_source")]
if (nrow(bad) > 0) {
  message("Different connectivity:")
  print(bad, row.names = FALSE, right = FALSE)
}
write.csv(bad, file.path(out, "check_rows_different_connectivity.csv"), row.names = FALSE)

fil <- rows[rows$role == "fill", ]
message("\nRows the change would fill:")
print(table(set = fil$set, status = fil$status))
message("Names sent in the shorter form: ", sum(fil$shortened), " of the fill rows, ",
        sum(fil$shortened & fil$status == "one_key"), " with one key.")


# ---- 2. What happens to the R01 pairs ----------------------------------------

# The same expert join as scripts/12_v0.7/02_acceptance.R. Empty labels are
# dropped, and a row without a label counts as not covered.
src <- "data-raw/comets"
norm_id <- function(x) toupper(gsub("[[:space:]]", "", x))
ids <- unlist(lapply(com$info$files$label, function(s)
  read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)[[1]]))
rk <- paste(com$rows$file, norm_id(ids), sep = "|")
mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character", na.strings = c("", "NA"))
mm <- mm[!is.na(mm$Manual_Metabolite_Group_Label), ]
mk <- paste(sub("\\.xlsx$", "", mm$SOURCE_FILE), norm_id(mm$IPT_METABID), sep = "|")
groups_of <- split(mm$Manual_Metabolite_Group_Label, mk)[rk]
names(groups_of) <- com$rows$row_id
covered <- function(a, b) lengths(groups_of[a]) > 0 & lengths(groups_of[b]) > 0
same <- function(a, b) vapply(seq_along(a), function(i)
  length(intersect(groups_of[[a[i]]], groups_of[[b[i]]])) > 0, logical(1))

simulate <- function(audit) {
  r <- audit$rows
  name_key <- stats::setNames(rows$name_key, rows$row_id)
  extra <- unname(name_key[r$row_id])
  key <- ifelse(r$ik_valid, r$inchikey, extra)
  from_name <- !r$ik_valid & !is.na(extra)
  p <- audit$pairs[audit$pairs$rule_id == "R01", c("row_a", "row_b", "group_id")]
  a <- match(p$row_a, r$row_id)
  b <- match(p$row_b, r$row_id)
  p$name_a <- r$original_name[a]
  p$name_b <- r$original_name[b]
  ka <- key[a]; kb <- key[b]
  same_name <- !is.na(r$name_normalized[a]) & r$name_normalized[a] == r$name_normalized[b]
  p$key_a <- ka; p$key_b <- kb
  p$outcome <- ifelse(is.na(ka) | is.na(kb), "no_key_on_a_side",
               ifelse(from_name[a] & from_name[b], "both_from_names",
               ifelse(ka == kb, "exact",
               ifelse(block1(ka) == block1(kb), "same_connectivity",
               ifelse(same_name, "different_connectivity_same_name", "different_connectivity")))))
  p
}

message("\n2. R01 pairs on the development data:")
pd <- simulate(dev)
print(pd[, c("name_a", "name_b", "outcome")], row.names = FALSE, right = FALSE)

pc <- simulate(com)
pc$covered <- covered(pc$row_a, pc$row_b)
pc$expert_same <- ifelse(pc$covered, same(pc$row_a, pc$row_b), NA)
message("\nR01 pairs on COMETS: ", nrow(pc))
print(table(outcome = pc$outcome, expert_same = pc$expert_same, useNA = "ifany"))

ex <- pc[pc$outcome == "exact", ]
ex_known <- ex[ex$covered, ]
share_exact <- if (nrow(ex_known) > 0) mean(ex_known$expert_same) else NA_real_
message("\n3. Would become exact: ", nrow(ex), " of ", nrow(pc), ". Needed: at least 80.")
message("4. Agree with the expert: ", sum(ex_known$expert_same), " of ", nrow(ex_known),
        " (", round(100 * share_exact, 1), "%). Needed: at least 98%.")
cols <- c("name_a", "name_b", "key_a", "key_b", "outcome", "expert_same")
if (any(!ex_known$expert_same)) {
  message("Exact against the expert:")
  print(ex_known[!ex_known$expert_same, cols], row.names = FALSE, right = FALSE)
}

di <- pc[pc$outcome == "different_connectivity" & pc$covered, ]
share_distinct <- if (nrow(di) > 0) mean(!di$expert_same) else NA_real_
allow_distinct <- nrow(di) > 0 && share_distinct >= 0.98
message("\nDifferent connectivity, different names: ", sum(pc$outcome == "different_connectivity"),
        "; covered by the expert: ", nrow(di), "; different compounds: ", sum(!di$expert_same),
        " (", round(100 * share_distinct, 1), "%). Distinct needs at least 98%: ",
        if (allow_distinct) "allowed" else "not allowed, they stay unresolved")
show_di <- pc[pc$outcome == "different_connectivity", ]
if (nrow(show_di) > 0) print(show_di[, cols], row.names = FALSE, right = FALSE)
write.csv(pc, file.path(out, "comets_r01_pairs.csv"), row.names = FALSE)
write.csv(pd, file.path(out, "development_r01_pairs.csv"), row.names = FALSE)


# ---- Decision ----------------------------------------------------------------

summary <- data.frame(
  condition = c("one key for rows with a key", "same connectivity", "COMETS R01 pairs becoming exact",
                "those pairs agree with the expert"),
  value = c(round(share_answered, 4), round(share_block1, 4), nrow(ex), round(share_exact, 4)),
  needed = c(0.5, 0.98, 80, 0.98))
summary$pass <- summary$value >= summary$needed
message("\n==== summary")
print(summary, row.names = FALSE, right = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
message("\nGo-ahead: ", if (isTRUE(all(summary$pass)) && n_err == 0) "YES" else "NO",
        if (n_err > 0) " (errors left, run again)" else "")
