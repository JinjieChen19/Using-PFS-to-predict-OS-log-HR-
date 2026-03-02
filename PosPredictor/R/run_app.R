#' Launch the PosPredictor Shiny Application
#'
#' Opens the four-tab Shiny application for computing the posterior
#' Probability of Success (PoS) of the OS endpoint using PFS as a surrogate.
#'
#' @param ... Arguments passed to \code{\link[shiny]{runApp}}.
#' @return Called for its side-effect of launching the Shiny app.
#' @export
run_app <- function(...) {
  app_dir <- system.file("shiny", package = "PosPredictor")
  if (!nchar(app_dir)) {
    stop("Could not find Shiny app directory. ",
         "Try re-installing 'PosPredictor'.")
  }
  shiny::runApp(app_dir, ...)
}
