
# Purpose: Perform DBSCAN clustering analysis on CNV data and analyze Disease-related gene regions not in clusters.
# === Helper Functions ===

# # Function to perform DBSCAN clustering on a dataframe using given eps and minPts
perform_dbscan <- function(df, eps, min_pts) {
  try(dbscan(df, eps = eps, minPts = min_pts, borderPoints =F))
}


# # Function to plot the DBSCAN clustering results and optionally highlight Disease-related gene regions
save_cluster_plot <- function(dat_pt, df_wt, cluster_centers, focal_calls_0_genes, focal_calls_0_dist_genes, index, plot_path, show_plot = T, save_plot=F) {
  suppressWarnings({
    p <- dat_pt %>%
      filter(DBscan_cluster != 0) %>%
      ggplot(aes(diagnosis_noise, relapse_noise)) +
      geom_abline(slope = 1, linetype = 2) +
      geom_point(data = df_wt, colour = "grey", alpha = 0.05, shape = 1) +
      geom_point(aes(colour = DBscan_cluster), alpha = 0.2) +
      geom_point(data = cluster_centers, colour = "red", size = 2) +
      geom_vline(xintercept = c(seq(0.2, 5, 1), seq(0.8, 4, 1)), linetype = 3) +
      geom_hline(yintercept = c(seq(0.2, 5, 1), seq(0.8, 4, 1)), linetype = 3) +
      {
        if (!is.null(focal_calls_0_genes) && length(focal_calls_0_genes) > 0) {
          geom_point(data = dat_pt %>% filter(gene_names %in% focal_calls_0_genes), colour = "blue", size = 2, aes(text = gene_names))
        }
      } +
      {
        if (!is.null(focal_calls_0_dist_genes) && length(focal_calls_0_dist_genes) > 0) {
          geom_point(data = dat_pt %>% filter(gene_names %in% focal_calls_0_dist_genes), colour = "green", size = 2, aes(text = gene_names))
        }
      } +
      ggrepel::geom_label_repel(data = dat_pt %>% filter(gene_names %in% focal_calls_0_dist_genes), aes(label = gene_names)) +
      geom_rect(aes(xmin = 1.8, xmax = 2.2, ymin = 1.8, ymax = 2.2), fill = "white", alpha = 0, color = "black") +
      coord_fixed() +
      ggtitle(paste0(index, " Diagnosis & Relapse")) +
      theme_light() +
      xlim(0, 4) + ylim(0, 4) +
      ylab("CN Relapse") + xlab("CN Diagnosis")
    
    if(save_plot){
    ggsave(plot_path, p, width = 10, height = 9, units = "in", dpi = 150)
    }
    
    if (show_plot) {
      print(p)
    }
    return(p)
    dev.off()
  })
}




# === Main DBSCAN Clustering & Focal Region Analysis Function ===
EvoCluster <- function(df, focalGR, plot_dir, pt_name, show_plot=F, save_plot=F) { 

  noise_sd <- 0.04  # standard deviation of noise added to points
  
  message("Processing patient: ", pt_name)
  
  # Identify column indices for diagnosis and relapse
  idx_d <- match(paste0(pt_name, "_D"), names(df))
  idx_r <- match(paste0(pt_name, "_R"), names(df))
  
  if (is.na(idx_d) | is.na(idx_r)) {
    stop("Columns for patient ", pt_name, " not found in dataset.")
  }

  # Prepare patient-specific dataframe with noise added to CN values
  dat_pt <- df %>%
    select(chromosome, start, end, gene_names, diagnosis = idx_d, relapse = idx_r) %>%
    mutate(diagnosis_noise = diagnosis + rnorm(nrow(df), mean = 0, sd = noise_sd),
           relapse_noise = relapse + rnorm(nrow(df), mean = 0, sd = noise_sd)) 
  setDT(dat_pt)
  
  
  # Identify data points in the "wild-type" (WT) region (CN ~ 2)
  df_wt <- dat_pt %>%
    filter(diagnosis_noise > 1.8 & diagnosis_noise < 2.2 &
             relapse_noise > 1.8 & relapse_noise < 2.2)
  setDT(df_wt)
  
  # Filter out the WT region from clustering
  dat_pt <- dat_pt %>%
    filter(!(diagnosis_noise > 1.8 & diagnosis_noise < 2.2 &
               relapse_noise > 1.8 & relapse_noise < 2.2))
  
  # Break data into blocks for DBSCAN (for performance)
  block_size <- 10000
  num_blocks <- ceiling(nrow(dat_pt) / block_size)
  
  dbs_clusters <- list()
  
  # Run DBSCAN on each block
  for (j in seq_len(num_blocks)) {
    
    start_row <- (j - 1) * block_size + 1
    end_row <- min(j * block_size, nrow(dat_pt))
    
    block_data <- dat_pt[start_row:end_row, ]
    
    dbs <- perform_dbscan(block_data %>% select(diagnosis_noise, relapse_noise), eps = 0.03, min_pts = 30)
    
    block_data$DBscan_cluster <- as.character(dbs$cluster)
    
    dbs_clusters[[j]] <- block_data

    rm(block_data, dbs)
    gc()
  }
  
  # Recombine all blocks after clustering
  dat_pt <- bind_rows(dbs_clusters)
  
  # Compute centers for each cluster
  cluster_centers <- dat_pt %>%
    filter(DBscan_cluster != "0") %>%
    group_by(DBscan_cluster) %>%
    summarise(diagnosis_noise = mean(diagnosis_noise),
              relapse_noise = mean(relapse_noise),
              n = n(),
              .groups = "drop")
  
  # Summarize number of genes per chromosome per cluster
  cluster_chromosomes <- dat_pt %>%
    group_by(DBscan_cluster, chromosome) %>%
    summarise(diagnosis_noise = mean(diagnosis_noise, na.rm = TRUE),
              relapse_noise = mean(relapse_noise, na.rm = TRUE),
              n_genes = n(),
              .groups = "drop") %>%
    mutate(pt_name = pt_name)
  
  # Check for overlaps between "noise points" and known disease-related gene regions
  if (!is.null(focalGR) & nrow(dat_pt) !=0) {
    dat_ptGR <- dat_pt  %>% filter(DBscan_cluster == 0) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # Keep only chromosomes present in both patient data and disease-related gene set
    chroms_in_focal <- unique(seqlevels(focalGR))
    chroms_in_focal <- intersect(as.character(seqlevels(dat_ptGR)), as.character(chroms_in_focal))
    
    dat_ptGR <- keepSeqlevels(dat_ptGR, chroms_in_focal, pruning.mode = "coarse")
    
    # Check for overlaps
    overlap_counts <- countOverlaps(focalGR, dat_ptGR)
    if (sum(overlap_counts) == 0) {
      message("No overlaps with Disease-Realted gene table in ", pt_name)
      focal_calls_0 <- tibble()
        focal_calls_0_dist <- tibble()
        focal_calls_0_dist_genes <- character(0)
    } else {
      overlaps <- findOverlaps(focalGR, dat_ptGR)
      focal_data_matched <- focalGR[queryHits(overlaps)] %>% as.data.frame()
      dat_data_matched <- dat_ptGR[subjectHits(overlaps)] %>% as.data.frame()
      focal_calls_0 <- bind_cols(dat_data_matched %>% select(-gene_names), select(focal_data_matched, CNA, gene_names))
    }
    
    # For matched focal genes, compute distance from all cluster centers
    if (nrow(focal_calls_0) > 0) {
      df1 <- focal_calls_0 %>% select(diagnosis_noise, relapse_noise)
      df2 <- cluster_centers %>% select(diagnosis_noise, relapse_noise)
      
      focal_calls_0_dist <- tryCatch({
        distmat <- dist(df1, df2, method = "euclidean")
        focal_calls_0$min_dist_clust_center <- apply(distmat, 1, min)
        focal_calls_0 %>% filter(min_dist_clust_center > 0.3)
      }, error = function(e) {
        message("Error computing distance matrix: ", e$message)
        tibble()
      })
      
      focal_calls_0_dist <- focal_calls_0_dist %>% mutate(pt_name = pt_name)
      focal_calls_0_dist_genes <- if (nrow(focal_calls_0_dist) > 0) focal_calls_0_dist$gene_names else character(0)
    } else {
      
      focal_calls_0_dist <- tibble()
      focal_calls_0_dist_genes <- character(0)
    }
  } else {
    focal_calls_0 <- tibble()
    focal_calls_0_dist <- tibble()
    focal_calls_0_dist_genes <- character(0)
  }
  
  # Plot clustering results
  p <- save_cluster_plot(
    dat_pt,
    df_wt,
    cluster_centers,
    if (exists("focal_calls_0") && nrow(focal_calls_0) > 0) focal_calls_0$gene_names else character(0),
    focal_calls_0_dist_genes,
    index = pt_name,
    file.path(plot_dir, paste0(pt_name, "_plot.png")),
    show_plot,
    save_plot
  )
  rm(dat_pt, df_wt)
  invisible(gc())
  
  message("Clustering complete - ", pt_name)
  
  # Return outputs
  list(
    all_clusters = cluster_chromosomes,
    all_focals_alone = focal_calls_0_dist,
    plot = p
  )
}

