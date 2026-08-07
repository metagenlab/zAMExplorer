#' Resolve an installed zAMPExplorer resource
#'
#' @param ... Path components inside the installed package.
#' @return A path to an installed package resource.
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "zAMPExplorer")
}
