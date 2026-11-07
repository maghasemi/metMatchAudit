# v0.11 feasibility (notes/decisions.md, 2026-11-07, v0.11 plan, step 2).
# Changes nothing in the package and sends nothing to PubChem. Works out
# which pairs the new rule R01s would take, on the v0.9.1 audits.
#
# R01s: an R01 pair whose two rows are in different files, have the same
# name after normalize_name(), and neither of whose names is a group name,
# a list of alternatives, an internal standard, an isomer label or marker,
# or a tentative name.
#
# Run from the project root:
#   Rscript scripts/18_v0.11/01_feasibility.R 2>&1 | tee results/v11_feasibility_log.txt

pkgload::load_all(".", quiet = TRUE)
stopifnot(as.character(utils::packageVersion("metMatchAudit")) == "0.9.1")
options(width = 200)

out <- "results/v11_feasibility"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

dev <- readRDS("results/v091_check/audit_development.rds")
com <- readRDS("results/v091_check/audit_comets.rds")
stopifnot(dev$info$package_version == "0.9.1", com$info$package_version == "0.9.1")


# ---- The rule's conditions on every pair -------------------------------------

flagged <- function(r) {
  cols <- c("name_group", "name_alternatives", "name_internal_standard", "name_isomer_label",
            "name_isomer_marker", "name_tentative")
  Reduce(`|`, lapply(cols, function(col) r[[col]] %in% TRUE))
}

pair_table <- function(audit, set) {
  r <- audit$rows
  p <- audit$pairs
  a <- match(p$row_a, r$row_id)
  b <- match(p$row_b, r$row_id)
  na <- r$name_normalized[a]
  nb <- r$name_normalized[b]
  data.frame(
    set = set, row_a = p$row_a, row_b = p$row_b, group_id = p$group_id,
    rule_id = p$rule_id, relationship = p$relationship,
    name_a = r$original_name[a], name_b = r$original_name[b],
    key_a = ifelse(r$ik_valid[a], r$inchikey[a], NA), key_b = ifelse(r$ik_valid[b], r$inchikey[b], NA),
    cross_file = r$file[a] != r$file[b],
    same_name = !is.na(na) & !is.na(nb) & na == nb,
    flagged = flagged(r)[a] | flagged(r)[b],
    keys_both = r$ik_valid[a] & r$ik_valid[b],
    stringsAsFactors = FALSE)
}
pt <- rbind(pair_table(dev, "development"), pair_table(com, "comets"))
pt$name_ok <- pt$cross_file & pt$same_name & !pt$flagged
pt$r01s <- pt$rule_id == "R01" & pt$name_ok


# ---- The expert groups on COMETS ---------------------------------------------

# The same join as scripts/12_v0.7/02_acceptance.R. Empty labels are
# dropped, and a row without a label counts as not covered.
src <- "data-raw/comets"
norm_id <- function(x) toupper(gsub("[[:space:]]", "", x))
ids <- unlist(lapply(com$info$files$label, function(s)
  read.csv(file.path(src, paste0(s, ".csv")), colClasses = "character", check.names = FALSE)[[1]]))
stopifnot(length(ids) == nrow(com$rows))
rk <- paste(com$rows$file, norm_id(ids), sep = "|")
mm <- read.csv(file.path(src, "ManualMappings.csv"), colClasses = "character", na.strings = c("", "NA"))
mm <- mm[!is.na(mm$Manual_Metabolite_Group_Label), ]
mk <- paste(sub("\\.xlsx$", "", mm$SOURCE_FILE), norm_id(mm$IPT_METABID), sep = "|")
groups_of <- split(mm$Manual_Metabolite_Group_Label, mk)[rk]
names(groups_of) <- com$rows$row_id

is_com <- pt$set == "comets"
pt$covered <- NA
pt$expert_same <- NA
pt$covered[is_com] <- lengths(groups_of[pt$row_a[is_com]]) > 0 & lengths(groups_of[pt$row_b[is_com]]) > 0
pt$expert_same[is_com] <- vapply(which(is_com), function(i)
  length(intersect(groups_of[[pt$row_a[i]]], groups_of[[pt$row_b[i]]])) > 0, logical(1))
pt$expert_same[is_com & !pt$covered] <- NA


# ---- 1 and 2. R01 pairs that would get same_name -------------------------------

r01 <- pt[pt$rule_id == "R01", ]
message("R01 pairs and the rule's conditions:")
print(table(set = r01$set,
            condition = ifelse(r01$r01s, "R01s",
                        ifelse(!r01$cross_file, "same file",
                        ifelse(!r01$same_name, "different names", "same name, flagged name")))))

message("\nSame name within one file (stay R01):")
show_cols <- c("set", "name_a", "name_b", "expert_same")
sf <- r01[!r01$cross_file & r01$same_name, ]
if (nrow(sf) > 0) print(sf[, show_cols], row.names = FALSE, right = FALSE) else message("  none")

message("\nSame name across files, but a flagged name (stay R01):")
fl <- r01[r01$cross_file & r01$same_name & r01$flagged, ]
if (nrow(fl) > 0) print(fl[, show_cols], row.names = FALSE, right = FALSE) else message("  none")

s <- pt[pt$r01s, ]
s_com <- s[s$set == "comets", ]
s_known <- s_com[s_com$covered, ]
share_expert <- if (nrow(s_known) > 0) mean(s_known$expert_same) else NA_real_
message("\nR01s pairs: development ", sum(s$set == "development"), ", COMETS ", nrow(s_com))
message("R01s pairs with a key on one side: ", sum(!is.na(s$key_a) | !is.na(s$key_b)), " of ", nrow(s))
message("1. COMETS R01 pairs that would get same_name: ", nrow(s_com), ". Needed: at least 80.")
message("2. In one expert group: ", sum(s_known$expert_same), " of ", nrow(s_known), " the expert covers (",
        round(100 * share_expert, 1), "%). Needed: at least 98%. Not covered: ", sum(!s_com$covered))
if (any(!s_known$expert_same)) {
  message("R01s pairs the expert separates:")
  print(s_known[!s_known$expert_same, c("name_a", "name_b", "key_a", "key_b", "group_id")],
        row.names = FALSE, right = FALSE)
}
if (sum(s$set == "development") > 0) {
  message("R01s pairs on the development data:")
  print(s[s$set == "development", c("name_a", "name_b", "key_a", "key_b")], row.names = FALSE, right = FALSE)
}
write.csv(s, file.path(out, "r01s_pairs.csv"), row.names = FALSE)


# ---- 3. Pairs with two keys that pass the name conditions --------------------

k <- pt[pt$keys_both & pt$name_ok, ]
k$conflict <- k$relationship == "distinct" | k$rule_id == "R12n"
share_ok <- if (nrow(k) > 0) mean(!k$conflict) else NA_real_
message("\n3. Pairs with two keys, different files and the same name: ", nrow(k),
        " (development ", sum(k$set == "development"), ", COMETS ", sum(k$set == "comets"), ")")
print(table(set = k$set, rule = k$rule_id))
message("Keys that do not conflict: ", sum(!k$conflict), " of ", nrow(k), " (", round(100 * share_ok, 1),
        "%). Needed: at least 98%.")
if (any(k$conflict)) {
  message("Pairs whose keys conflict:")
  print(k[k$conflict, c("set", "name_a", "name_b", "key_a", "key_b", "rule_id", "expert_same")],
        row.names = FALSE, right = FALSE)
}
k_com <- k[k$set == "comets" & k$covered, ]
if (nrow(k_com) > 0) {
  message("For comparison, the same pairs on COMETS against the expert: ", sum(k_com$expert_same), " of ",
          nrow(k_com), " in one group")
}
write.csv(k, file.path(out, "keyed_same_name_pairs.csv"), row.names = FALSE)


# ---- Decision ----------------------------------------------------------------

summary <- data.frame(
  condition = c("COMETS R01 pairs getting same_name", "those pairs in one expert group",
                "keyed same-name pairs without a key conflict"),
  value = c(nrow(s_com), round(share_expert, 4), round(share_ok, 4)),
  needed = c(80, 0.98, 0.98))
summary$pass <- summary$value >= summary$needed
message("\n==== summary")
print(summary, row.names = FALSE, right = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
message("\nGo-ahead: ", if (isTRUE(all(summary$pass))) "YES" else "NO")
