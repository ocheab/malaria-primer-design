#' @import shiny
#' @import shinyWidgets
#' @import bslib
#' @import Biostrings
#' @import DT
#' @import ggplot2
#' @import seqinr
#' @import stringr
NULL
# In R/run_app.R
run_app <- function() {
  shiny::shinyApp(ui = ui, server = server)
}
