#' Launch the Malaria Primer Designer App
#' @export
run_app <- function() {
  addResourcePath("www", system.file("www", package = "malariaPrimerDesigner"))
  shiny::shinyApp(ui = ui, server = server)
}
