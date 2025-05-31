library(shiny)
library(Biostrings)
library(DT)
library(stringr)
library(ggplot2)
library(seqinr)
library(shinyWidgets)
library(bslib)
library(shinybusy)

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

ui <- tagList(
  add_busy_spinner(spin = "fading-circle", color = "#2C3E50", position = "bottom-right"),
  navbarPage(
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
))

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
  fasta_set <- if (!is.null(input$customFasta)) {
    readDNAStringSet(input$customFasta$datapath, format = "fasta")
  } else {
    req(input$marker, input$species)
    file_path <- file.path(root_dir, input$marker, paste0(input$species, ".txt"))
    readDNAStringSet(file_path, format = "fasta")
  }
  
  req(fasta_set)
  
  primers <- data.frame(
    Accession = character(), Species = character(), Gene = character(),
    Forward = character(), Reverse = character(), Start = integer(), End = integer(),
    Amplicon_Size = integer(),
    Fwd_Tm = numeric(), Rev_Tm = numeric(), Delta_Tm = numeric(),
    Fwd_GC = numeric(), Rev_GC = numeric(), GC = numeric(),
    Hairpin = character(), Dimer = character(), GC_Clamp = character(),
    Score = numeric(), stringsAsFactors = FALSE
  )
  
  withProgress(message = "Designing primers...", value = 0, {
    total_seqs <- length(fasta_set)
    max_pairs_per_seq <- 10
    rev_window <- c(100, 1000)  # search up to 600 bp downstream
    
    for (j in seq_along(fasta_set)) {
      incProgress(j / total_seqs, detail = paste("Processing sequence", j, "of", total_seqs))
      
      # Only extract first token of FASTA header
      this_acc <- strsplit(names(fasta_set)[j], "\\s+")[[1]][1]
      this_seq <- gsub("\\s", "", as.character(fasta_set[[j]]))
      seqlen <- nchar(this_seq)
      found_pairs <- 0
      
      for (i in 1:(seqlen - input$maxLen - rev_window[2])) {
        for (fwd_len in input$minLen:input$maxLen) {
          if ((i + fwd_len - 1) > seqlen) next
          fwd <- substr(this_seq, i, i + fwd_len - 1)
          fwd_gc <- (str_count(fwd, "[GC]") / fwd_len) * 100
          fwd_tm <- calculate_tm(fwd)
          
          if (fwd_gc < input$gcRange[1] || fwd_gc > input$gcRange[2] ||
              fwd_tm < input$tmRange[1] || fwd_tm > input$tmRange[2]) next
          
          # Sample reverse primer region every 10 bp
          for (k in seq(i + rev_window[1], min(i + rev_window[2], seqlen - input$minLen), by = 10)) {
            for (rev_len in input$minLen:input$maxLen) {
              if ((k + rev_len - 1) > seqlen) break
              rev_seq <- substr(this_seq, k, k + rev_len - 1)
              rev <- as.character(reverseComplement(DNAString(rev_seq)))
              rev_gc <- (str_count(rev, "[GC]") / rev_len) * 100
              rev_tm <- calculate_tm(rev)
              
              if (rev_gc < input$gcRange[1] || rev_gc > input$gcRange[2] ||
                  rev_tm < input$tmRange[1] || rev_tm > input$tmRange[2]) next
              
              amplicon_size <- (k + rev_len - 1) - i + 1
              
              # Relaxed hairpin rule: allow <= 6 matches with self
              relaxed_hairpin_check <- function(seq) {
                rc <- as.character(reverseComplement(DNAString(seq)))
                score <- sum(strsplit(seq, "")[[1]] == strsplit(rc, "")[[1]])
                return(score > 6)
              }
              
              hp <- ifelse(relaxed_hairpin_check(fwd) | relaxed_hairpin_check(rev), "Yes", "No")
              dm <- ifelse(has_cross_dimer(fwd, rev), "Yes", "No")
              clamp <- ifelse(has_gc_clamp(fwd) & has_gc_clamp(rev), "Yes", "No")
              
              score <- abs(fwd_tm - mean(input$tmRange)) +
                abs(rev_tm - mean(input$tmRange)) +
                abs(fwd_gc - mean(input$gcRange)) +
                abs(rev_gc - mean(input$gcRange)) +
                ifelse(hp == "Yes", 2, 0) +
                ifelse(dm == "Yes", 2, 0)
              
              primers <- rbind(primers, data.frame(
  Accession = this_acc,
  Species = ifelse(!is.null(input$customFasta), "User FASTA", input$species),
  Gene = ifelse(!is.null(input$customFasta), "Custom", input$marker),
  Forward = fwd, Reverse = rev,
  Start = i, End = k + rev_len - 1,
  Amplicon_Size = amplicon_size,
  Fwd_Tm = round(fwd_tm, 2), Rev_Tm = round(rev_tm, 2),
  Delta_Tm = round(abs(fwd_tm - rev_tm), 2),
  Fwd_GC = round(fwd_gc, 2), Rev_GC = round(rev_gc, 2),
  GC = round((fwd_gc + rev_gc) / 2, 2),
  Hairpin = hp, Dimer = dm, GC_Clamp = clamp,
  Score = round(score, 2),
  Mode = "Strict"  # or "Relaxed"
))
              
              found_pairs <- found_pairs + 1
              if (found_pairs >= max_pairs_per_seq) break
            }
            if (found_pairs >= max_pairs_per_seq) break
          }
          if (found_pairs >= max_pairs_per_seq) break
        }
        if (found_pairs >= max_pairs_per_seq) break
      }
    }
    
    if (nrow(primers) == 0) {
  showNotification("No primers found with strict constraints. Relaxing parameters...", type = "warning")

  for (j in seq_along(fasta_set)) {
    incProgress(j / total_seqs, detail = paste("Relaxed search for sequence", j))
    
    this_acc <- strsplit(names(fasta_set)[j], "\\s+")[[1]][1]
    this_seq <- gsub("\\s", "", as.character(fasta_set[[j]]))
    seqlen <- nchar(this_seq)
    found_pairs <- 0
    
    for (i in 1:(seqlen - input$maxLen - rev_window[2])) {
      for (fwd_len in input$minLen:input$maxLen) {
        if ((i + fwd_len - 1) > seqlen) next
        fwd <- substr(this_seq, i, i + fwd_len - 1)
        fwd_gc <- (str_count(fwd, "[GC]") / fwd_len) * 100
        fwd_tm <- calculate_tm(fwd)

        # Relax: no GC or Tm filtering
        for (k in seq(i + rev_window[1], min(i + rev_window[2], seqlen - input$minLen), by = 10)) {
          for (rev_len in input$minLen:input$maxLen) {
            if ((k + rev_len - 1) > seqlen) break
            rev_seq <- substr(this_seq, k, k + rev_len - 1)
            rev <- as.character(reverseComplement(DNAString(rev_seq)))
            rev_gc <- (str_count(rev, "[GC]") / rev_len) * 100
            rev_tm <- calculate_tm(rev)
            
            amplicon_size <- (k + rev_len - 1) - i + 1
            
            hp <- "Unknown"
            dm <- "Unknown"
            clamp <- "Unknown"
            
            score <- abs(fwd_tm - mean(input$tmRange)) +
              abs(rev_tm - mean(input$tmRange)) +
              abs(fwd_gc - mean(input$gcRange)) +
              abs(rev_gc - mean(input$gcRange))

            primers <- rbind(primers, data.frame(
              Accession = this_acc,
              Species = ifelse(!is.null(input$customFasta), "User FASTA", input$species),
              Gene = ifelse(!is.null(input$customFasta), "Custom", input$marker),
              Forward = fwd, Reverse = rev,
              Start = i, End = k + rev_len - 1,
              Amplicon_Size = amplicon_size,
              Fwd_Tm = round(fwd_tm, 2), Rev_Tm = round(rev_tm, 2),
              Delta_Tm = round(abs(fwd_tm - rev_tm), 2),
              Fwd_GC = round(fwd_gc, 2), Rev_GC = round(rev_gc, 2),
              GC = round((fwd_gc + rev_gc) / 2, 2),
              Hairpin = hp, Dimer = dm, GC_Clamp = clamp,
              Score = round(score, 2)
            ))

            found_pairs <- found_pairs + 1
            if (found_pairs >= max_pairs_per_seq) break
          }
          if (found_pairs >= max_pairs_per_seq) break
        }
        if (found_pairs >= max_pairs_per_seq) break
      }
      if (found_pairs >= max_pairs_per_seq) break
    }
  }

  showNotification(paste("Relaxed search completed:", nrow(primers), "primers found."), type = "message")
}    
    
    primers <- primers[order(primers$Score), ]
    primer_results(primers)
    output$primerTable <- renderDT({
  df <- primer_results()
  req(df)
  df$Mode <- ifelse(df$Mode == "Relaxed",
                  '<span style="color: white; background-color: #e67e22; padding: 2px 6px; border-radius: 4px;">Relaxed</span>',
                  '<span style="color: white; background-color: #2ecc71; padding: 2px 6px; border-radius: 4px;">Strict</span>')
  datatable(df, selection = 'multiple', filter = 'top',
            options = list(pageLength = 10), rownames = FALSE)
})
    
    showNotification(paste("Primer design completed:", nrow(primers), "primers found."), type = "message")
    
    output$primerTable <- renderDT({
      df <- primer_results()
      req(df)
      datatable(df, selection = 'multiple', filter = 'top',
                options = list(pageLength = 10), rownames = FALSE)
    })
  })
})
  
  output$tmPlot <- renderPlot({
  shinybusy::show_modal_spinner(text = "Generating Tm distribution plot...", spin = "fading-circle")
  
  req(primer_results())
  df <- primer_results()
  
  # Convert wide Tm to long format
  df_long <- tidyr::pivot_longer(df, cols = c(Fwd_Tm, Rev_Tm), names_to = "Type", values_to = "Tm")
  
  p <- ggplot(df_long, aes(x = Tm, fill = Type)) +
    geom_histogram(position = "dodge", binwidth = 1, color = "black") +
    theme_minimal() +
    labs(title = "Tm Distribution", x = "Melting Temperature (°C)", y = "Count") +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    scale_fill_manual(values = c("Fwd_Tm" = "#F8766D", "Rev_Tm" = "#00BFC4"))
  
  shinybusy::remove_modal_spinner()
  p
})
  
  output$gcPlot <- renderPlot({
  shinybusy::show_modal_spinner(text = "Generating GC content plot...", spin = "fading-circle")
  
  req(primer_results())
  df <- primer_results()

  # ✅ Make sure the columns exist and are numeric
  if (!("Fwd_GC" %in% names(df)) || !("Rev_GC" %in% names(df))) {
    shinybusy::remove_modal_spinner()
    return(NULL)
  }

  df_long_gc <- tidyr::pivot_longer(
    df,
    cols = c("Fwd_GC", "Rev_GC"),
    names_to = "Type",
    values_to = "GC_Value"
  )

  # ✅ Check for valid rows
  if (nrow(df_long_gc) == 0 || all(is.na(df_long_gc$GC_Value))) {
    shinybusy::remove_modal_spinner()
    showNotification("No valid GC content to plot.", type = "error")
    return(NULL)
  }

  p <- ggplot(df_long_gc, aes(x = GC_Value, fill = Type)) +
    geom_histogram(position = "dodge", binwidth = 1, color = "black") +
    theme_minimal() +
    labs(title = "GC Content Distribution", x = "GC%", y = "Count") +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    scale_fill_manual(values = c("Fwd_GC" = "#7CAE00", "Rev_GC" = "#C77CFF"))

  shinybusy::remove_modal_spinner()
  p
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
