#' metMatchAudit: audit cross-dataset metabolite matches
#'
#' Matching tools such as MetLinkR put rows from different datasets under
#' one standardized name. A shared name does not always mean a shared
#' structure: it can hide positional isomers, stereoisomers, a broader
#' annotation next to a specific one, or two unrelated compounds. This
#' package looks at the structural evidence for each pair of rows in a
#' group and reports how they relate, together with the rule that
#' decided it.
#'
#' The rules and their order are described in `notes/design.md` in the
#' project repository.
#'
#' @keywords internal
"_PACKAGE"
