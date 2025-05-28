#' Launch the Malaria Primer Designer Shiny App
#'
#' This function launches the interactive user interface for designing PCR primers
#' for malaria diagnostic markers across Plasmodium species.
#'
#' @export
run_app <- function() {
  shiny::shinyApp(ui = ui, server = server)
}
