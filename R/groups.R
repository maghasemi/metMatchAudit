# Most severe first. A group is summarised by its most severe pair.
# same_name compared no structure, so it ranks below every structural
# answer except unresolved.
.severity <- c("distinct", "unresolved", "same_name", "broader_narrower", "stereo_variant", "exact")

#' Summarise pair results per group
#'
#' @param pairs Output of [classify_pairs()].
#'
#' @return A data frame with one line per group: `group_id`, `n_rows`,
#'   `n_pairs`, `worst_relationship` (most severe first: distinct,
#'   unresolved, same_name, broader_narrower, stereo_variant, exact), `worst_rule`
#'   (the rule behind the first pair with that relationship), and one
#'   count column per relationship.
#'
#' @export
summarise_groups <- function(pairs) {
  needed <- c("group_id", "row_a", "row_b", "relationship", "rule_id")
  missing <- setdiff(needed, names(pairs))
  if (length(missing) > 0) {
    stop("`pairs` is missing ", paste(missing, collapse = ", "),
         ". Run classify_pairs() first.", call. = FALSE)
  }
  unknown <- setdiff(unique(pairs$relationship), .severity)
  if (length(unknown) > 0) {
    stop("Unknown relationship in `pairs`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }

  groups <- unique(pairs$group_id)
  rows <- lapply(groups, function(g) {
    p <- pairs[pairs$group_id == g, , drop = FALSE]
    rank <- match(p$relationship, .severity)
    worst <- which(rank == min(rank))[1]
    counts <- vapply(.severity, function(r) sum(p$relationship == r), integer(1))
    data.frame(group_id = g,
               n_rows = length(unique(c(p$row_a, p$row_b))),
               n_pairs = nrow(p),
               worst_relationship = p$relationship[worst],
               worst_rule = p$rule_id[worst],
               t(counts),
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  out <- if (length(rows) > 0) do.call(rbind, rows) else {
    data.frame(group_id = character(), n_rows = integer(), n_pairs = integer(),
               worst_relationship = character(), worst_rule = character(),
               stringsAsFactors = FALSE)
  }
  names(out)[names(out) %in% .severity] <- paste0("n_", names(out)[names(out) %in% .severity])
  rownames(out) <- NULL
  out
}
