#' Launch zAMPExplorer
#'
#' Launch the zAMPExplorer Shiny application for downstream exploration of
#' 16S amplicon microbiome data stored in a `phyloseq` object.
#'
#' @param onStart Optional function called before the application starts.
#' @param options Named list of Shiny application options.
#' @param enableBookmarking Bookmarking mode passed to [shiny::shinyApp()].
#' @param uiPattern URL pattern passed to [shiny::shinyApp()].
#' @param ... Reserved for backward compatibility.
#'
#' @return A `shiny.appobj`.
#' @export
#' @examples
#' if (interactive()) {
#'   zAMPExplorer_app()
#' }
zAMPExplorer_app <- function(
    onStart = NULL,
    options = list(),
    enableBookmarking = NULL,
    uiPattern = "/",
    ...) {
  shiny::shinyApp(
    ui = app_ui,
    server = app_server,
    onStart = onStart,
    options = options,
    enableBookmarking = enableBookmarking,
    uiPattern = uiPattern
  )
}
