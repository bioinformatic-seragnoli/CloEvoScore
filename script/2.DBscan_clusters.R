
# Purpose: Perform DBSCAN clustering analysis on CNV data and analyze Disease-related gene regions not in clusters.
# === Helper Functions ===

# Perform DBSCAN clustering
perform_dbscan <- function(df, eps, min_pts) {
  dbscan(df, eps = eps, minPts = min_pts, borderPoints = FALSE)
}


# Trajectory Plot 
save_cluster_plot <- function(dat_pt, df_wt, cluster_centers, focal_calls_0_genes, focal_calls_0_dist_genes, index, plot_path, show_plot = T) {
  suppressWarnings({
    p <- dat_pt %>%
      filter(DBscan_cluster != 0) %>%
      ggplot(aes(diagnosis_noise, relapse_noise)) +
      geom_abline(slope = 1, linetype = 2) +
      ggforce::geom_mark_hull(aes(group = DBscan_cluster), expand = 0.008) +
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
      ggtitle(paste0("Patient n°", index, " Diagnosis & Relapse")) +
      theme_light() +
      xlim(0, 4) + ylim(0, 4) +
      ylab("CN Relapse") + xlab("CN Diagnosis")
    
    
    ggsave(plot_path, p, width = 10, height = 9, units = "in", dpi = 300)
    
    if (show_plot) {
      print(p)
    }
    return(p)
  })
}





# === Main Analysis ===
EvoCluster <- function(df,eps_pair, min_pair, noise_sd, focalGR, plot_dir, pt_list, show_plot=F) {  
  set.seed(123)
  all_clusters <- list()
  all_focals_alone <- list()
  plot_list <- list()
  
  # Loop through each patient in the list
  for (i in seq_along(pt_list)) {
    message("Procesing patient index: ", i)
    
    idx_d <- (i * 2 - 1) + 4
    idx_r <- (i * 2) + 4
    
    pt_name <- pt_list[i]
    
    # Prepare patient-specific dataset
    dat_pt <- df %>%
      select(chromosome, start, end, gene_names, diagnosis = idx_d, relapse = idx_r) %>%
      mutate(diagnosis_noise = diagnosis + rnorm(nrow(df), mean = 0, sd = noise_sd),
             relapse_noise = relapse + rnorm(nrow(df), mean = 0, sd = noise_sd)) 
    
    df_wt <-
      dat_pt %>%  filter((
        diagnosis_noise > 1.8 &
          diagnosis_noise < 2.2 & relapse_noise > 1.8 & relapse_noise < 2.2
      ))
    
    
    dat_pt <-
      dat_pt %>%  filter(!(
        diagnosis_noise > 1.8 &
          diagnosis_noise < 2.2 & relapse_noise > 1.8 & relapse_noise < 2.2
      ))
    
    # Perform DBSCAN clustering
    dbs <- perform_dbscan(dat_pt %>% select(diagnosis_noise, relapse_noise), eps = eps_pair, min_pts = min_pair)
    dat_pt$DBscan_cluster <- as.character(dbs$cluster)
    
    # Summarize cluster centers
    cluster_centers <- dat_pt %>%
      filter(DBscan_cluster != 0) %>%
      group_by(DBscan_cluster) %>%
      summarise(diagnosis_noise = mean(diagnosis_noise),
                relapse_noise = mean(relapse_noise),
                n = n(),
                .groups = "drop")
    
    # Counting the number of genes contributing to each cluster for each chromosome
    cluster_chromosomes <- dat_pt %>%
      group_by(DBscan_cluster, chromosome) %>%
      summarise(
        diagnosis_noise = mean(diagnosis_noise, na.rm = TRUE),
        relapse_noise = mean(relapse_noise, na.rm = TRUE),
        n_genes = n(),
        .groups = "drop"
      )
    
    cluster_chromosomes$pt_idx <- i
    cluster_chromosomes$pt_name <- pt_name
    
    cluster_chromosomes$pt_idx <- i
    cluster_chromosomes$pt_name <- pt_name
    all_clusters[[i]] <- cluster_chromosomes
    
    # If Disease-related regions are provided, perform overlap analysis
    # Perform overlap analysis if focal regions are provided
    focal_calls_0 <- tibble()
    if (!is.null(focalGR)) {
      
    # Match focal regions
    dat_ptGR <- makeGRangesFromDataFrame(dat_pt, keep.extra.columns = T)
    
    suppressWarnings(seqlevelsStyle(focalGR) <- "NCBI")
    suppressWarnings(seqlevelsStyle(dat_ptGR) <- "NCBI")
    
    #Maintaining only common levels
    common_levels <- intersect(seqlevels(focalGR), seqlevels(dat_ptGR))
    focalGR <- keepSeqlevels(focalGR, common_levels, pruning.mode = "coarse")
    dat_ptGR <- keepSeqlevels(dat_ptGR, common_levels, pruning.mode = "coarse")
    
    focal_calls <- plyranges::join_overlap_intersect(focalGR, dat_ptGR) %>% as.data.frame()

    focal_calls_0 <- focal_calls %>% filter(DBscan_cluster == 0)
    
    df1 = focal_calls_0 %>% select(diagnosis_noise, relapse_noise)
    df2 = cluster_centers %>% select(diagnosis_noise, relapse_noise)
    
    # Creating a distance matrix from Disease-related genes and cluster
    try({
      distmat <- dist(df1, df2, method="euclidean")
      focal_calls_0$min_dist_clust_center <- apply(distmat, 1, min) 
      focal_calls_0_dist <- focal_calls_0 %>% filter(min_dist_clust_center> 0.3)
      })
    
    
    tryCatch({
      focal_calls_0_dist$pt_idx <- i
      focal_calls_0_dist$pt_name <- pt_name
      
      focal_calls_0$pt_idx <- i
      focal_calls_0$pt_name <- pt_name
    }, error = function(e) {
      message("The patient does not have Disease-related genes outside clusters.")
      focal_calls_0_dist <- tibble()  # Return an empty tibble if there is an error
    })
    
    
    if (exists("focal_calls_0") && nrow(focal_calls_0) > 0) {
      focal_calls_0_dist_genes <- focal_calls_0_dist$gene_names
    } else {
      focal_calls_0_dist_genes <- character(0)
    }
    } else {
      focal_calls_0_dist_genes <- character(0)
      focal_calls_0_dist <- tibble()
    }
    
    all_focals_alone[[i]] <- focal_calls_0_dist
    
    # Making Trajectory plot 
    p <- save_cluster_plot(
      dat_pt,
      df_wt,
      cluster_centers,
      focal_calls_0$hgnc_symbol,
      focal_calls_0_dist_genes,
      index = i,
      file.path(plot_dir, paste0(pt_name, "_plot.png")),
      show_plot
      )
    plot_list[[pt_name]] <- p
    }
    
    
    message(paste("Clustering complete - pt", pt_name))
  

  list(all_clusters = bind_rows(all_clusters),
       all_focals_alone =  bind_rows(all_focals_alone),
       plot_list = plot_list)
}



