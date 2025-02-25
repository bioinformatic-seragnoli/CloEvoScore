# Load required libraries
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(furrr) 
library(dbscan)


server <- function(input, output, session) {
  
  observe({
    showModal(
      modalDialog(
        title =  tags$div(
          style = "background-color: #1b4f72; color: #ffffff; padding: 15px; margin: -15px -15px 15px -15px;text-align: center;",
          tagList(
          tags$img(src = "logo2_1.png", height = "50px", style ="vertical-align: middle; margin-right: 10px;"),
          tags$span("Welcome to CloEvoScore app", style = "vertical-align: middle;")
        )),
        HTML("
          <p>This interactive open source tool calculates the evolutionary trajectories of patients based on gene copy number values.</p>
          
          <p>The app enables rapid and intuitive analysis of patients' genomic evolution by allowing you to:</p>
          <ul>
            <li>Upload gene copy number data.</li>
            <li>Compute evolutionary trajectories and derive evolutionary scores.</li>
            <li>Visualize results through interactive plots and summary tables.</li>
            <li>Download the results for further analysis.</li>
          </ul>
          
          <p>For more details on the pipeline, please refer to the reference paper: 
          <a href='https://doi.org/10.1101/347534' target='_blank'>Preprint available here</a>.</p>
          
          <p>If you need any help with input files, please click the <strong>?</strong> icon in the header.</p>
          
          <p><em>IMPORTANT: This example uses genomic profiles of 5 multiple myeloma samples created solely for demonstrating the online version of this tool.</em></p>
        "),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  observeEvent(input$help, {
    showModal(modalDialog(
      title = div(
        style = "background-color: #1b4f72; color: #ffffff; padding: 15px; margin: -15px -15px 15px -15px; text-align: center;",
        tagList(icon("file-code"), "Input File Requirements")
      ),
      p("This application requires specific input files for accurate analysis. Below are the descriptions and requirements for each file:"),
      
      tags$ul(
        tags$li(
          tags$b("Gene Call Table"),
          " contains detailed information on the copy number values of individual genes for each patient in the study. These values represent the number of gene copies present in the genome and are essential for understanding genetic alterations associated with the disease being analyzed. This file serves as the foundation for identifying patterns of genomic amplification or deletion across the patient cohort.",
          tags$ul(
            tags$li(tags$b("chromosome:"), " The chromosome where the gene is located."),
            tags$li(tags$b("start:"), " Start position of the gene."),
            tags$li(tags$b("end:"), " End position of the gene."),
            tags$li(tags$b("gene_names:"), " Gene symbol."),
            tags$li("At least two columns per sample with copy number values at two time points (e.g., 'PAZ_001_D', 'PAZ_001_R').")
          )
        ),
        tags$li(
          tags$b("Ploidy and Purity Table"),
          " (Optional) provides two critical parameters used to refine the interpretation of copy number values. ",
          "Ploidy refers to which value to center the normality region of the copy number profile, while purity indicates the proportion of cancer cells within a given sample compared to normal cells. Correcting the raw copy number values using ploidy and purity ensures a more accurate reflection of true genetic alterations.",
          tags$ul(
            tags$li(tags$b("pt_name:"), " Patient’s code (must match the Gene Call Table)."),
            tags$li(tags$b("Ploidy_D"), " and ", tags$b("Ploidy_R:"), " Ploidy values for each time point."),
            tags$li(tags$b("Purity_D"), " and ", tags$b("Purity_R:"), " Purity values for each time point.")
          )
        ),
        tags$li(
          tags$b("Disease-Related Gene Regions"),
          " (Optional) contains a curated list of genes that are particularly relevant to the disease under investigation. These genes are typically identified through computational analyses using tools such as ",
          tags$b("GISTIC2.0"),
          ". In the analysis, if these genes are not associated with any specific cluster, they will still be considered individually and contribute to the final score.",
          tags$ul(
            tags$li(tags$b("gene_name:"), " Gene symbol."),
            tags$li(tags$b("chr:"), " Chromosome where the gene is located."),
            tags$li(tags$b("start:"), " Start position of the gene."),
            tags$li(tags$b("end:"), " End position of the gene."),
            tags$li(tags$b("CNA:"), " Copy number alteration detected by GISTIC (e.g., 'A' for Amplification, 'D' for Deletion).")
          )
        ),
        tags$li("All files must be provided in ", tags$b("CSV, TSV, or TXT"), " format.")
      ),
      
      p("Examples about file structure can be found on the application's Github page. If you need further assistance, please contact us at bioinformatic.seragnoli@gmail.com."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  
  
  
  noise_sd <- 0.04
  eps_pair <- 0.03
  min_pair <- 30
  set.seed(123)
  
  validate_input <- function(file, required_cols) {
    data <- fread(file$datapath)
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      stop(paste("Input file is missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    return(data)
  }
  
  load_data <- function(input_file, default_path, required_cols = NULL) {
    if (!is.null(input_file)) {
      cat("Loading user-provided file...\n")
      data <- fread(input_file$datapath)
      if (!is.null(required_cols)) {
        missing_cols <- setdiff(required_cols, names(data))
        if (length(missing_cols) > 0) {
          stop(paste("File is missing required columns:", paste(missing_cols, collapse = ", ")))
        }
      }
      return(data)
    } else {
      cat("Loading default file from", default_path, "...\n")
      return(fread(file.path("www", default_path)))
    }
  }
  
  
  
  cat("Setting up directories for outputs and plots...\n")
  output_dir <- "outputs"
  plot_dir <- file.path(output_dir, "plots/")
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(plot_dir, showWarnings = FALSE)
  
  analysisComplete <- reactiveVal(FALSE)
  finalScores <- reactiveVal(NULL)
  plotList <- reactiveVal(NULL)
  
  observeEvent(input$run_analysis, {
    
    analysisComplete(FALSE) 
    output$status <- renderText("Analysis in progress...")
    shinyjs::hide("patient_select")
    
    
    cat("Reading input files or using defaults...\n")
    gene_calls <- load_data(input$gene_calls, "data/example/gene_table.txt", c("chromosome", "start", "end", "gene_names"))
    pp_table <- if (input$toggle_ploidy_purity) {
      load_data(input$ploidy_purity, "data/example/test_pp.txt", c("pt_name", "Purity_D", "Purity_R", "Ploidy_D", "Ploidy_R"))
    } else NULL
    focal_data <- if (input$toggle_disease_regions) {
      load_data(input$disease_regions, "data/focal_loci_hg19.txt", c("chromosome", "start", "end", "gene_name"))
    } else NULL
    
    
    pt_list <- names(gene_calls)[-c(1:4)] %>% 
      str_remove_all(pattern = paste0(c("_[Rr]$", "_[Er]$", "_[Dd]$"), collapse = "|")) %>% unique()
    
    cat("Loading functions and EVO-score tables...\n")
    pos_evol <- fread("../data/POS_pressures.txt")
    neg_evol <- fread("../data/NEG_pressures.txt")
    
    source("../script/1.Purity_and_ploidy_gene_correction.R")
    source("../script/2.DBscan_clusters.R")
    source("../script/3.EVOscore_assignment.R")
    source("../script/4.Trajectory_estimation.R")
    
    progressr::with_progress({
      p <- progressr::progressor(steps = 4)
      
      if (!is.null(pp_table)) {
        p("Step 1: Correcting CN Values")
        cat("Step 1: Correcting CN Values...\n")
        corrected_calls <- correct_cn(gene_calls, pp_table)
        write_tsv(corrected_calls, file.path(output_dir, "1.geneCN_CALLS_samples_FINAL.txt"))
      } else {
        corrected_calls <- gene_calls
      }
      
      p("Step 2: Performing Clustering")
      cat("Step 2: Performing Clustering...\n")
      if (!is.null(focal_data)) {
        focalGR <- makeGRangesFromDataFrame(focal_data, keep.extra.columns = TRUE)
        clustering <- EvoCluster(corrected_calls, eps_pair, min_pair, noise_sd, focalGR, plot_dir, pt_list, show_plot=F)
        write_tsv(clustering$all_clusters, file.path(output_dir, "2.dbscan_results.txt"))
        write_tsv(clustering$all_focals_alone, file.path(output_dir, "2.focals_not_in_clusters_results.txt"))
      } else {
        clustering <- EvoCluster(corrected_calls, eps_pair, min_pair, noise_sd, focalGR =NULL, plot_dir, pt_list, show_plot=F) 
        write_tsv(clustering$all_clusters, file.path(output_dir, "2.dbscan_results.txt"))
      }
      
      
      p("Step 3: Assigning Evolutionary Scores")
      cat("Step 3: Assigning Evolutionary Scores...\n")
      state_thresholds <- tibble(
        states = c("hT", "T", "AT", "A", "sA", "wt", "sD", "D", "HD", "H"),
        limits_low = c(4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2, -Inf),
        limits_high = c(Inf, 4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2)
      )
      cluster_results <- assign_labels(clustering$all_clusters, state_thresholds)
      evo_clusters <- evoC_scores(cluster_results, pos_evol, neg_evol)
      write_tsv(evo_clusters, file.path(output_dir, "3.cluster_evoscores.txt"))
      
      if (!is.null(focal_data)) {
        focal_results <- assign_labels(clustering$all_focals_alone, state_thresholds)
        evo_gene <- evoF_scores(focal_results, pos_evol, neg_evol)
        write_tsv(evo_gene, file.path(output_dir, "3.focal_evoscores.txt"))
      } else {
        evo_gene <- tibble(pt_name=pt_list)
      }
      
      
      p("Step 4: Calculating Trajectories")
      cat("Step 4: Calculating Trajectories...\n")
      final_scores <- CloEvoScore(evo_clusters, evo_gene, output_dir)
      write_tsv(final_scores, file.path(output_dir, "4.Trajectories.txt"))
    })
    
    cat("Analysis complete!\n")
    analysisComplete(TRUE)
    patient_names <-  pt_list
    updateSelectInput(session, "patient_select", choices = patient_names)
    output$status <- renderText("Analysis complete! Select Patient to visualize his trajectory")
    
    finalScores(final_scores) 
    plotList(clustering$plot_list)
    shinyjs::show("patient_select")
  })
  
  
  output$analysisComplete <- renderUI({
    if (analysisComplete()) {
      tags$div("Analysis is complete!")
    } else {
      NULL
    }
  })
  
  output$results_box <- renderUI({
    div(id = "resbox",
      box(
    title =  tagList(icon("chart-line"), "Results"),
    status = "success",
    solidHeader = F,
    width = 12,
    collapsible = T,
    collapsed  = !analysisComplete(),
    class = "centered-box",
    textOutput("status"),
    selectInput("patient_select", "Select Patient Results to Visualize:", choices = NULL),
    fluidRow(
      box(
        title = "Trajectory Plot",
        height = "600px",
        width = 8,
        collapsible = F,
        solidHeader = TRUE,
        status = "success",
        class = "plot-container",
        plotOutput("progress_plot", height = "550px", width = "100%")
      ),
      box(
        title = "EVO-scores and trajectory",
        height = "400px",
        width = 4,
        collapsible = F,
        solidHeader = TRUE,
        status = "primary",
        tableOutput("summary_results")
      )
    ),
    
  ))})
  
  output$download_button <- renderUI({
    if (analysisComplete()) {
      fluidRow(
        column(
          width = 12,
          downloadButton("download_results", "Download Results")
        )
      )
    } else {
      NULL
    }
  })
  
  
  
  observeEvent(input$patient_select, {
    
    output$progress_plot <- renderPlot({
      req(input$patient_select)
      selected_patient <- input$patient_select
      plot_list <- plotList()
      
      if (!is.null(plot_list) && selected_patient %in% names(plot_list)) {
        p <- plot_list[[selected_patient]] 
        par(mar=c(4,4,2,1))
        print(p)
      } else {
        
        plot(1, type = "n", xlab = "", ylab = "", axes = FALSE)
        text(1, 1, "No plot available for the selected patient.", cex = 1.5)
      }
      
    })
    
    
    output$summary_results <- renderTable({
      req(finalScores()) 
      req(input$patient_select)
      
      selected_patient <- input$patient_select
      filtered_data <- finalScores() %>% filter(pt_name == selected_patient)
      
      
      transposed_data <- as.data.frame(t(filtered_data))
      transposed_data <-rownames_to_column(transposed_data, "Info")
      
      colnames(transposed_data) <- c("Info", "Value")
      transposed_data
    })
    
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste("final_results-", Sys.Date(), ".zip", sep = "")  # Dynamic filename with date
    },
    content = function(file) {
      cat("Preparing download...\n")
      zip::zipr(file, files = list.files(output_dir, full.names = TRUE))  # Create a zip archive of outputs
    }
  )
  

  
}
