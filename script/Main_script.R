#!/usr/bin/env Rscript

# === Environment Setup ======================================
# Load required libraries
library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(purrr) 
library(dbscan)
library(factoextra)
library(proxy)
library(ggrepel)
library(ggforce)

message("Starting script execution...")

# === Define Command Line Options =============================
option_list <- list(
  make_option(c("-g", "--genes"), type = "character", default = NULL,
              help = "Gene copy number table [default: %default]"),
  make_option(c("-p", "--pp"), type = "character", default = NULL,
              help = "Purity and ploidy table [default: %default]"),
  make_option(c("-r", "--rel"), type = "character", default = NULL,
              help = "Disease-related genes [default: %default]"),
  
  make_option(c("-o", "--output"), type = "character", default = "outputs/",
              help = "Output directory [default: %default]"),
  make_option(c("-f", "--plots"), type = "character", default = "plots/",
              help = "Plots directory [default: %default]")
)


# Parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))
message("Parsed command line arguments.")

# === Configure Directories ============================
data_dir <- "/cloevoscore/data" #please change directory path in 'www/data' if you execute the pipeline outside the docker
output_dir <- file.path(opt$output)
plot_dir <- file.path(output_dir, opt$plots)
message("Output directory set to: ", output_dir)
message("Plots directory set to: ", plot_dir)

# Create output directories if they do not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
message("Directories checked and created if necessary.")

# === Load External Scripts ===================================
script_dir <- "script/"
message("Loading external scripts...")
source(file.path(script_dir, "load_files.R"))
source(file.path(script_dir, "1.Purity_and_ploidy_gene_correction.R"))
source(file.path(script_dir, "2.DBscan_clusters.R"))
source(file.path(script_dir, "3.EVOscore_assignment.R"))
source(file.path(script_dir, "4.Trajectory_estimation.R"))
message("External scripts loaded successfully.")


# === Verify Input Files ======================================
if (!file.exists(opt$genes)) stop("Error: Gene copy number file does not exist.")
if (!file.exists(opt$pp)) stop("Error: Purity and ploidy file does not exist.")
if (!file.exists(opt$rel)) stop("Error: Disease-related genes file does not exist.")
message("All input files verified.")

# === Load Input Data ========================================
message("Loading input data...")
data <- load_data(gene_tb = opt$genes, pp_tb = opt$pp, rel_genes = opt$rel)
pt_list <- data$pt_list
message("Data loaded successfully. Patients identified: ", length(pt_list))


res <- list()
for (pt in pt_list) {
  
  message(paste("Processing patient", pt))
  
  sample_dir <- paste0(output_dir,"/", pt)
  dir.create(sample_dir)
  
  # Identify relevant columns for the patient
  column_names <- names(data$gene_calls)
  pt_cols <- column_names[grepl(paste0("^", pt, "_"), column_names)]
  selected_cols <- c("chromosome", "start", "end", "gene_names", pt_cols)
  
  # Extract patient-specific gene calls
  gene_calls <- data$gene_calls %>%  select(all_of(selected_cols)) %>% as.data.table()
  
# Apply purity and ploidy correction if data is available
if(!is.null(data$pp)){
  message("Applying purity and ploidy correction for patient: ", pt)
corrected_calls <- correct_cn(gene_calls, data$pp %>% filter(pt_name == pt))
write.csv(corrected_calls, paste0(sample_dir, "/genecalls_ploidypurity_corr.txt"))
} else {
  corrected_calls <- gene_calls
  rm(gene_calls)
  gc()
}

  # Perform clustering analysis if focal data is available
if(!is.null(data$rel_genes)){
  message("Running clustering analysis for patient: ", pt)
focalGR <- makeGRangesFromDataFrame(data$rel_genes, keep.extra.columns = T)

clustering_results <- EvoCluster(corrected_calls, focalGR, plot_dir, pt, show_plot=F, save_plot=T)
if (nrow(clustering_results$all_focals_alone) !=0) {
  write.csv(clustering_results$all_focals_alone, paste0(sample_dir, "/diseaserelated_genes_off.txt"))
}
} else {
  clustering_results <- EvoCluster(corrected_calls, focalGR=NULL, plot_dir, pt, show_plot=F, save_plot=T)
}
  
  message("Clustering complete for patient: ", pt)
  rm(corrected_calls)
  gc()
  write.csv(clustering_results$all_clusters, paste0(sample_dir, "/dbscan_results.txt"))
 
  
# Assign cluster labels
message(paste("Assigning scores for:", pt))
state_thresholds <- tibble(
  states =  c("hT","T","AT","A","sA","wt","sD","D","HD","H"),
  limits_low = c(4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2, -Inf),
  limits_high = c(Inf, 4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2)
)


# Process clusters
cluster_labels <- assign_labels(clustering_results$all_clusters, state_thresholds)
evo_clusters <- evoC_scores(cluster_labels, data$pos_evol_tab, data$neg_evol_tab)

if (!is.null(data$rel)) {
  focal_labels <- assign_labels(as.data.frame(clustering_results$all_focals_alone), state_thresholds)
  evo_gene <- evoF_scores(focal_labels, data$pos_evol_tab, data$neg_evol_tab)
} else {
  evo_gene <- tibble(pt_name=pt)
}

message("Calculating trajectory for patient: ", pt)
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

message(paste("Process complete for patient:", pt))

}

# Save final scores
all_scores <- bind_rows(lapply(res, function(x) if (!is.null(x$final_scores)) x$final_scores else NULL), .id = "patient")
write.csv(all_scores, paste0(output_dir, "/trajectories_and_scores.txt"))

message("Analysis complete! Results saved to ", output_dir)

