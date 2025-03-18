source('packages.R')
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
          <a href='https://i.imgflip.com/qz3zz.jpg' target='_blank'>Preprint available here</a>.</p>
          
          <p>If you need any help with input files, please click the <strong>?</strong> icon in the header.</p>
          
          <p><em>IMPORTANT: This example uses genomic profiles of 4 multiple myeloma samples created solely for demonstrating the online version of this tool.</em></p>
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
        tags$img(src = "logo2_1.png", height = "50px", style ="vertical-align: middle; margin-right: 10px;"),
        tags$span("Input File Requirements", style = "vertical-align: middle;")
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
            tags$li("At least two columns per sample with copy number values at two time points (e.g., 'PT_001_D', 'PT_001_R').")
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
  
  # --------------------------- Data prep --------------------------------------
  source("script/shinyapp_ausiliar_functions.R")
  cat("Setting up directories for outputs and plots...\n")
  output_dir <- "outputs"
  plot_dir <- file.path(output_dir, "plots/")
  dir.create(output_dir, showWarnings = FALSE)
  
  
  analysisComplete <- reactiveVal(FALSE)
  finalScores <- reactiveVal(NULL)
  plotList <- reactiveVal(NULL)
  
  observeEvent(input$run_analysis, {
    
    send_status_update <- function(session, message) {
      session$sendCustomMessage(type = "statusUpdate", message)
    }
    
    analysisComplete(FALSE) 
    send_status_update(session, "Analysis in progress...")
    shinyjs::hide("patient_select")
    
    withProgress(message = "Running analysis...", value = 0, {
      
      
      cat("Reading input files or using defaults...\n")
      pp_table <- if (input$toggle_ploidy_purity) {
        load_data(input$ploidy_purity, "data/example/test_pp.txt", c("pt_name", "Purity_D", "Purity_R", "Ploidy_D", "Ploidy_R"))
      } else NULL
      focal_data <- if (input$toggle_disease_regions) {
        load_data(input$disease_regions, "data/focal_loci_hg19.txt", c("chromosome", "start", "end", "gene_names"))
      } else NULL
      

      path_gene_calls <- ifelse(is.null(input$gene_calls), "www/data/example/gene_table.txt", input$gene_calls$datapath)
      column_names <- names(fread(path_gene_calls, nrows = 1, header = TRUE, data.table = FALSE))
      
      pt_list <- column_names[!column_names %in% c("chromosome", "start", "end", "gene_names")] %>% 
        str_remove_all(pattern = paste0("_[RrEeDd]$", collapse = "|")) %>% unique()
      
      # ---------------------------------------- CloEvoScore pipeline --------------------------------
      
      cat("Loading functions and EVO-score tables...\n")
      
      pos_evol <- fread("www/data/POS_pressures.txt")
      neg_evol <- fread("www/data/NEG_pressures.txt")
      
      source("script/1.Purity_and_ploidy_gene_correction.R")
      source("script/2.DBscan_clusters.R")
      source("script/3.EVOscore_assignment.R")
      source("script/4.Trajectory_estimation.R")
      
      res <- list()
      for (pt in pt_list) {
        sample_dir <- paste0(output_dir,"/", pt)
        dir.create(sample_dir)
        pt_cols <- column_names[grepl(paste0("^", pt, "_"), column_names)]

        selected_cols <- c("chromosome", "start", "end", "gene_names", pt_cols)
        gene_calls <- fread(path_gene_calls, select = selected_cols, data.table = FALSE)
        setDT(gene_calls)
        
        if (!is.null(pp_table) && pt %in% pp_table$pt_name) {
        incProgress(0.25 / length(pt_list), detail = paste("Correcting CN for:", pt))
        
        send_status_update(session, "Correcting CN Values...")
        p("Step 1: Correcting CN Values")
        cat("Step 1: Correcting CN Values...\n")
     
        corrected_calls <- correct_cn(gene_calls, pp_table %>% filter(pt_name == pt))
        write.csv(corrected_calls, paste0(sample_dir, "/genecalls_ploidypurity_corr.txt"))
        rm(gene_calls)
        gc()
      } else {
        corrected_calls <- gene_calls
        rm(gene_calls)
        gc()
        
      }
    
      
      p("Step 2: Performing Clustering")
      cat("Step 2: Performing Clustering...\n")
      set.seed(123)
      incProgress(0.25 / length(pt_list), detail = paste("Performing clustering for:", pt))
      if (!is.null(focal_data)) {
        focalGR <- makeGRangesFromDataFrame(focal_data, keep.extra.columns = TRUE)

       clustering_results <- EvoCluster(corrected_calls, focalGR, plot_dir, pt, show_plot=F, save_plot=F)
       if (nrow(clustering_results$all_focals_alone) !=0) {
         write.csv(clustering_results$all_focals_alone, paste0(sample_dir, "/diseaserelated_genes_off.txt"))
       }
        message("Clustering completed for: ", pt)
      } else {
  
       clustering_results <- EvoCluster(corrected_calls, focalGR=NULL, plot_dir, pt, show_plot=F, save_plot=F)
        message("Clustering completed for: ", pt)
       
      }
      rm(corrected_calls)
      gc()
      write.csv(clustering_results$all_clusters, paste0(sample_dir, "/dbscan_results.txt"))
      
      p("Step 3: Assigning Evolutionary Scores")
      cat("Step 3: Assigning Evolutionary Scores...\n")
      incProgress(0.25 / length(pt_list), detail = paste("Assigning scores for:", pt))
      state_thresholds <- tibble(
        states = c("hT", "T", "AT", "A", "sA", "wt", "sD", "D", "HD", "H"),
        limits_low = c(4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2, -Inf),
        limits_high = c(Inf, 4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2)
      )
      cluster_labels <- assign_labels(clustering_results$all_clusters, state_thresholds)
      evo_clusters <- evoC_scores(cluster_labels, pos_evol, neg_evol)
      
      if (!is.null(focal_data)) {
        focal_labels <- assign_labels(as.data.frame(clustering_results$all_focals_alone), state_thresholds)
        evo_gene <- evoF_scores(focal_labels, pos_evol, neg_evol)
      } else {
        evo_gene <- tibble(pt_name=pt)
      }
    
      
      p("Step 4: Calculating Trajectories")
      cat("Step 4: Calculating Trajectories...\n")
      final_scores <- CloEvoScore(evo_clusters, evo_gene, output_dir)
      
      
      res[[pt]] <- list(
        cluster_results = clustering_results,
        evo_clusters = evo_clusters,
        evo_gene = evo_gene,
        final_scores = final_scores,
        plot = if ("plot" %in% names(clustering_results)) clustering_results$plot else NULL
      )
      
      rm(clustering_results, evo_clusters, evo_gene, final_scores)
      gc()  
      
      incProgress(0.25 / length(pt_list), detail = paste("Completed:", pt))
      }
      
      plot_list <- lapply(res, function(x) if (!is.null(x$plot)) x$plot else NULL)
      final_scores <- bind_rows(lapply(res, function(x) if (!is.null(x$final_scores)) x$final_scores else NULL), .id = "patient")
      write.csv(final_scores, paste0(output_dir, "/trajectories_and_scores.txt"))
      
      send_status_update(session, "Analysis complete! Select Patient to visualize his trajectory")
      cat("Analysis complete!\n")  
      
      analysisComplete(TRUE)
      patient_names <-  pt_list
      updateSelectInput(session, "patient_select", choices = patient_names)
      output$status <- renderText("Analysis complete! Select Patient to visualize his trajectory")
      
      finalScores(final_scores) 
      plotList(plot_list)
      shinyjs::show("patient_select")
      
      
    })
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
              height = "450px",
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
      gc()
    }, res=96)
    
    
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
