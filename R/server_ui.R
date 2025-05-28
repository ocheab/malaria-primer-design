library(shiny)
library(Biostrings)
library(DT)
library(stringr)
library(ggplot2)
library(seqinr)
library(shinyWidgets)
library(bslib)

has_self_complementarity <- function(seq) {
  rc <- as.character(reverseComplement(DNAString(seq)))
  score <- sum(strsplit(seq, "")[[1]] == strsplit(rc, "")[[1]])
  return(score > 4)
}

has_gc_clamp <- function(seq) {
  tail <- substr(seq, nchar(seq) - 1, nchar(seq))
  return(grepl("[GC]{1,2}$", tail))
}

has_cross_dimer <- function(fwd, rev) {
  rev_rc <- as.character(reverseComplement(DNAString(rev)))
  for (len in 4:8) {
    fwd_tail <- substr(fwd, nchar(fwd) - len + 1, nchar(fwd))
    rev_head <- substr(rev_rc, 1, len)
    pairs <- mapply(function(a, b) {
      paste0(a, b) %in% c("AT", "TA", "GC", "CG")
    }, strsplit(fwd_tail, "")[[1]], strsplit(rev_head, "")[[1]])
    if (all(pairs)) return(TRUE)
  }
  return(FALSE)
}

calculate_tm <- function(seq) {
  seq <- toupper(seq)
  nA <- str_count(seq, "A")
  nT <- str_count(seq, "T")
  nG <- str_count(seq, "G")
  nC <- str_count(seq, "C")
  if (nchar(seq) < 14) {
    return(2 * (nA + nT) + 4 * (nG + nC))
  } else {
    return(64.9 + 41 * (nG + nC - 16.4) / nchar(seq))
  }
}

root_dir <- system.file("extdata", package = "malariaPrimerDesigner")
marker_choices <- basename(list.dirs(root_dir, recursive = FALSE))

ui <- navbarPage(
  title = div(
    img(src = "www/logo.png", height = "40px", style = "margin-right: 10px;"),
    span("Centre for Malaria and Other Tropical Diseases Care, UITH, Ilorin, Nigeria", style = "font-size: 16px; font-weight: bold;")
  ),
  theme = bs_theme(version = 5, bootswatch = "flatly", base_font = font_google("Lato")),
  
  tabPanel("Home",
           fluidRow(
             column(12,
                    h3("Welcome to the Malaria Primer Designer", style = "margin-top: 20px;"),
                    p("This tool helps you design specific and optimized primers for malaria diagnostic markers across Plasmodium species.", style = "font-size: 16px;"),
                    br(),
                    img(src = "www/pcr_diagram.png", height = "300px", style = "display: block; margin-left: auto; margin-right: auto;"),
                    p(em("Figure: Overview of the PCR process showing denaturation, primer annealing, and extension."), style = "text-align: center;"),
                    br(),
                    tags$ul(
                      tags$li(strong("Step 1 - Denaturation:"), " DNA is heated to separate strands."),
                      tags$li(strong("Step 2 - Primer Annealing:"), " Short primers bind to the target sequences."),
                      tags$li(strong("Step 3 - Extension:"), " DNA polymerase synthesizes new DNA strands from the primers.")
                    )
             )
           )
  ),
  
  tabPanel("Design Primers",
           sidebarLayout(
             sidebarPanel(
               fileInput("customFasta", "Upload Custom FASTA File (optional)", accept = ".fasta"),
               selectInput("marker", "Select Diagnostic Marker:", choices = marker_choices),
               uiOutput("speciesUI"),
               uiOutput("customModeNote"),
               numericInput("minLen", "Min Primer Length:", 18, min = 10),
               numericInput("maxLen", "Max Primer Length:", 25, min = 10),
               sliderInput("gcRange", "GC Content Range (%)", min = 30, max = 70, value = c(40, 60)),
               sliderInput("tmRange", "Melting Temperature (Tm °C)", min = 45, max = 70, value = c(55, 65)),
               actionButton("designBtn", label = "Design Primers", icon = icon("flask"), class = "btn-primary"),
               downloadButton("downloadPrimers", "Download Filtered Table as CSV"),
               downloadButton("downloadFasta", "Download Selected Primers as FASTA")
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Primer Table", DTOutput("primerTable", height = "500px")),
                 tabPanel("Tm Distribution", plotOutput("tmPlot")),
                 tabPanel("GC Content Distribution", plotOutput("gcPlot")),
                 tabPanel("Target Sequence", verbatimTextOutput("targetSeq"))
               )
             )
           )
  ),
  
  tabPanel("Help",
           fluidRow(
             column(12,
                    h4("Help & Instructions"),
                    tags$ul(
                      tags$li("Upload a FASTA file if you want to design primers on a custom sequence."),
                      tags$li("Otherwise, select a diagnostic marker and species with reference sequences available."),
                      tags$li("Adjust the GC content and Tm ranges to optimize primer quality."),
                      tags$li("Download your primers in table (CSV) or FASTA format for lab use.")
                    ),
                    p("For support, contact: cemtrod.ilorin@gmail.com, +2348140699446")
             )
           )
  )
)

server <- function(input, output, session) {
  addResourcePath("www", system.file("www", package = "malariaPrimerDesigner"))
  output$speciesUI <- renderUI({
    if (!is.null(input$customFasta)) {
      return(NULL)  # hide species selection when custom FASTA is uploaded
    }
    req(input$marker)
    marker_path <- file.path(root_dir, input$marker)
    species_files <- list.files(marker_path, pattern = "\\.txt$", full.names = FALSE)
    species_names <- gsub(".txt$", "", species_files)
    selectInput("species", "Select Plasmodium Species:", choices = species_names)
  })
  
  output$customModeNote <- renderUI({
    if (!is.null(input$customFasta)) {
      helpText("You are working with a custom FASTA upload. Built-in marker/species selection is disabled.")
    }
  })
  
  get_target_sequence <- reactive({
    if (!is.null(input$customFasta)) {
      fasta_set <- readDNAStringSet(input$customFasta$datapath, format = "fasta")
    } else {
      req(input$marker, input$species)
      file_path <- file.path(root_dir, input$marker, paste0(input$species, ".txt"))
      if (!file.exists(file_path)) return(NULL)
      fasta_set <- readDNAStringSet(file_path, format = "fasta")
    }
    paste(as.character(fasta_set), collapse = "")
  })
  
  output$targetSeq <- renderText({
    seq <- get_target_sequence()
    if (is.null(seq)) return("No sequence available.")
    return(seq)
  })
  
  primer_results <- reactiveVal()
  
  observeEvent(input$designBtn, {
    seq <- get_target_sequence()
    req(seq)
    seq <- gsub("\\s", "", seq)
    
    withProgress(message = "Designing primers...", value = 0, {
      incProgress(0.1, detail = "Initializing")
      seqlen <- nchar(seq)
      primers <- data.frame(
        Species = character(), Gene = character(), Forward = character(), Reverse = character(),
        Start = integer(), End = integer(), Length = integer(),
        GC = numeric(), Tm = numeric(), Hairpin = character(), Dimer = character(),
        GC_Clamp = character(), Score = numeric(), stringsAsFactors = FALSE
      )
      
      incProgress(0.3, detail = "Searching candidate primers")
      for (i in 1:(seqlen - input$maxLen)) {
        for (len in input$minLen:input$maxLen) {
          if ((i + len - 1) > seqlen) break
          
          fwd <- substr(seq, i, i + len - 1)
          rev <- as.character(reverseComplement(DNAString(fwd)))
          
          # GC and Tm for both forward and reverse
          fwd_gc <- (str_count(fwd, "[GC]") / len) * 100
          rev_gc <- (str_count(rev, "[GC]") / len) * 100
          fwd_tm <- calculate_tm(fwd)
          rev_tm <- calculate_tm(rev)
          
          # Check both primers fall within filter criteria
          if (fwd_gc >= input$gcRange[1] && fwd_gc <= input$gcRange[2] &&
              rev_gc >= input$gcRange[1] && rev_gc <= input$gcRange[2] &&
              fwd_tm >= input$tmRange[1] && fwd_tm <= input$tmRange[2] &&
              rev_tm >= input$tmRange[1] && rev_tm <= input$tmRange[2]) {
            
            hp <- ifelse(has_self_complementarity(fwd) | has_self_complementarity(rev), "Yes", "No")
            dm <- ifelse(has_cross_dimer(fwd, rev), "Yes", "No")
            clamp <- ifelse(has_gc_clamp(fwd) & has_gc_clamp(rev), "Yes", "No")
            
            score <- abs(fwd_tm - mean(input$tmRange)) +
              abs(rev_tm - mean(input$tmRange)) +
              abs(fwd_gc - mean(input$gcRange)) +
              abs(rev_gc - mean(input$gcRange)) +
              ifelse(hp == "Yes", 2, 0) +
              ifelse(dm == "Yes", 2, 0)
            
            primers <- rbind(primers, data.frame(
              Species = ifelse(!is.null(input$customFasta), "User FASTA", input$species),
              Gene = ifelse(!is.null(input$customFasta), "Custom", input$marker),
              Forward = fwd, Reverse = rev, Start = i, End = i + len - 1, Length = len,
              GC = round((fwd_gc + rev_gc) / 2, 2),  # Average GC
              Tm = round((fwd_tm + rev_tm) / 2, 2),  # Average Tm
              Hairpin = hp, Dimer = dm, GC_Clamp = clamp,
              Score = round(score, 2)
            ))
          }
        }
      }
      
      incProgress(0.9, detail = "Finalizing")
      if (nrow(primers) == 0) {
        showNotification("No primers found with current parameters. Try relaxing constraints.", type = "warning")
      }
      
      primers <- primers[order(primers$Score), ]
      primer_results(primers)
      showNotification(paste("Primer design completed:", nrow(primers), "primers found."), type = "message")
      output$primerTable <- renderDT({
        datatable(primers, selection = 'multiple', filter = 'top', options = list(pageLength = 10), rownames = FALSE)
      })
    })
  })
  
  output$tmPlot <- renderPlot({
    req(primer_results())
    ggplot(primer_results(), aes(x = Tm)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
      theme_minimal() + labs(title = "Melting Temperature Distribution", x = "Tm (°C)", y = "Count")
  })
  
  output$gcPlot <- renderPlot({
    req(primer_results())
    ggplot(primer_results(), aes(x = GC)) +
      geom_histogram(binwidth = 1, fill = "darkgreen", color = "black") +
      theme_minimal() + labs(title = "GC Content Distribution", x = "GC%", y = "Count")
  })
  
  output$downloadPrimers <- downloadHandler(
    filename = function() { paste0("filtered_primers_", Sys.Date(), ".csv") },
    content = function(file) {
      # Use the DataTable’s proxy to get filtered data
      dt_proxy <- isolate(input$primerTable_rows_all)  # Get indices of visible (filtered) rows
      all_primers <- primer_results()
      if (!is.null(dt_proxy)) {
        filtered_data <- all_primers[dt_proxy, ]
        write.csv(filtered_data, file, row.names = FALSE)
      } else {
        write.csv(all_primers, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadFasta <- downloadHandler(
    filename = function() {
      paste0("selected_primers_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      req(input$primerTable_rows_selected)
      primers <- primer_results()[input$primerTable_rows_selected, ]
      fasta_lines <- c()
      for (i in 1:nrow(primers)) {
        fasta_lines <- c(fasta_lines,
                         paste0(">", primers$Species[i], "_", primers$Gene[i], "_F", i),
                         primers$Forward[i],
                         paste0(">", primers$Species[i], "_", primers$Gene[i], "_R", i),
                         primers$Reverse[i])
      }
      writeLines(fasta_lines, con = file)
    }
  )
}

shinyApp(ui = ui, server = server)
