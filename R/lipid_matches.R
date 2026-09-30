#' Propose lipid matches across files
#'
#' Pairs every lipid row with every lipid row from a different file that
#' has the same lipid Species, applies the pair rules, and keeps the
#' pairs that come out `exact` or `broader_narrower`. It does not need a
#' matching step such as MetLinkR; if the rows carry a group column from
#' one, the result says which proposals that step had already made.
#'
#' A row can match several rows in the other file, for example the Broad
#' Species "PC 34:4" and the Metabolon structures "14:0/20:4" and
#' "16:1/18:3". Such proposals are marked `one_to_many`. The tool does not
#' pick one: a Species-level value usually covers all of them, which
#' matters when comparing measurements.
#'
#' @param evidence Rows after [row_evidence()] and [lipid_evidence()],
#'   with a `file` column.
#' @param group_col Optional name of a column with an existing grouping,
#'   such as MetLinkR's standardized name. `NULL` to skip.
#'
#' @return A data frame with one line per proposal: `row_a`, `row_b`,
#'   `file_a`, `file_b`, `name_a`, `name_b`, `lipid_species`,
#'   `relationship`, `rule_id`, `evidence`, `partners_a` (how many rows
#'   `row_a` is matched to), `partners_b`, `shape` (`one_to_one` or
#'   `one_to_many`) and, when `group_col` is given, `already_grouped`.
#'
#' @export
propose_lipid_matches <- function(evidence, group_col = "group_id") {
  needed <- c("row_id", "file", "original_name", "lipid_species")
  missing <- setdiff(needed, names(evidence))
  if (length(missing) > 0) {
    stop("`evidence` is missing ", paste(missing, collapse = ", "),
         ". Run row_evidence() and lipid_evidence() first.", call. = FALSE)
  }
  if (!is.null(group_col) && !group_col %in% names(evidence)) {
    stop("`evidence` has no column called `", group_col, "`.", call. = FALSE)
  }

  empty <- data.frame(row_a = character(), row_b = character(),
                      file_a = character(), file_b = character(),
                      name_a = character(), name_b = character(),
                      lipid_species = character(), relationship = character(),
                      rule_id = character(), evidence = character(),
                      partners_a = integer(), partners_b = integer(),
                      shape = character(), stringsAsFactors = FALSE)
  if (!is.null(group_col)) empty$already_grouped <- logical()

  lip <- evidence[!is.na(evidence$lipid_species), c("row_id", "file", "lipid_species")]
  files <- unique(evidence$file)
  pieces <- list()
  for (i in seq_along(files)) for (j in seq_along(files)) {
    if (j <= i) next
    a <- lip[lip$file == files[i], ]
    b <- lip[lip$file == files[j], ]
    m <- merge(a, b, by = "lipid_species", suffixes = c("_a", "_b"))
    if (nrow(m) > 0) pieces[[length(pieces) + 1]] <- m
  }
  if (length(pieces) == 0) return(empty)
  cand <- do.call(rbind, pieces)

  pairs <- data.frame(group_id = cand$lipid_species, row_a = cand$row_id_a,
                      row_b = cand$row_id_b, stringsAsFactors = FALSE)
  pairs <- classify_pairs(evidence, pairs)
  keep <- pairs$relationship %in% c("exact", "broader_narrower")
  if (!any(keep)) return(empty)
  p <- pairs[keep, , drop = FALSE]

  idx_a <- match(p$row_a, evidence$row_id)
  idx_b <- match(p$row_b, evidence$row_id)
  count_a <- table(p$row_a)
  count_b <- table(p$row_b)

  out <- data.frame(
    row_a = p$row_a, row_b = p$row_b,
    file_a = evidence$file[idx_a], file_b = evidence$file[idx_b],
    name_a = evidence$original_name[idx_a], name_b = evidence$original_name[idx_b],
    lipid_species = p$group_id,
    relationship = p$relationship, rule_id = p$rule_id, evidence = p$evidence,
    partners_a = as.integer(count_a[p$row_a]),
    partners_b = as.integer(count_b[p$row_b]),
    stringsAsFactors = FALSE
  )
  out$shape <- ifelse(out$partners_a == 1 & out$partners_b == 1, "one_to_one", "one_to_many")
  if (!is.null(group_col)) {
    ga <- evidence[[group_col]][idx_a]
    gb <- evidence[[group_col]][idx_b]
    out$already_grouped <- !is.na(ga) & !is.na(gb) & ga == gb
  }
  rownames(out) <- NULL
  out
}
