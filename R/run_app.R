# In R/run_app.R
run_app <- function() {
  shiny::shinyApp(ui = ui, server = server)
}