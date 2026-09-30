#' All pairs of rows within each group
#'
#' @param rows A data frame with `row_id` and a group column.
#' @param group_col Name of the column that holds the group, `"group_id"`
#'   by default. Rows with a missing group are left out.
#'
#' @return A data frame with columns `group_id`, `row_a` and `row_b`
#'   (row IDs). A group of n rows gives n * (n - 1) / 2 pairs, in input
#'   order. Groups with one row give none.
#'
#' @export
make_pairs <- function(rows, group_col = "group_id") {
  for (col in c("row_id", group_col)) {
    if (!col %in% names(rows)) {
      stop("`rows` needs a column called `", col, "`.", call. = FALSE)
    }
  }
  if (anyDuplicated(rows$row_id)) {
    stop("`row_id` must be unique.", call. = FALSE)
  }
  groups <- as.character(rows[[group_col]])
  keep <- !is.na(groups) & nzchar(groups)
  ids_by_group <- split(rows$row_id[keep], factor(groups[keep], levels = unique(groups[keep])))

  pieces <- lapply(names(ids_by_group), function(g) {
    ids <- ids_by_group[[g]]
    if (length(ids) < 2) return(NULL)
    idx <- utils::combn(length(ids), 2)
    data.frame(group_id = g, row_a = ids[idx[1, ]], row_b = ids[idx[2, ]],
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, pieces)
  if (is.null(out)) {
    out <- data.frame(group_id = character(), row_a = character(),
                      row_b = character(), stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  out
}
