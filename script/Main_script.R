#!/usr/bin/env Rscript

# === Environment Setup ======================================
# Load required libraries
library(optparse)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(furrr) 
library(dbscan)
library(factoextra)
library(proxy)
library(ggrepel)
library(ggforce)

# Define command line options
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


# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Configure paths
data_dir <- "data/"
script_dir <- "script/"
output_dir <- opt$output
plot_dir <- file.path(output_dir, opt$plots)

# Create output directories if they do not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load scripts
source(file.path(script_dir, "load_files.R"))
source(file.path(script_dir, "1.Purity_and_ploidy_gene_correction.R"))
source(file.path(script_dir, "2.DBscan_clusters.R"))
source(file.path(script_dir, "3.EVOscore_assignment.R"))
source(file.path(script_dir, "4.Trajectory_estimation.R"))


# === Constants ===
noise_sd <- 0.04
eps_pair <- 0.03
min_pair <- 30
set.seed(123)


data <- load_data(gene_tb = opt$genes, pp_tb = opt$purityploidy, rel_genes = opt$rel)

if(!is.null(data$pp)){
corrected_calls <- correct_cn(data$gene_calls, data$pp)
write_tsv(corrected_calls, file.path(output_dir, "1.geneCN_CALLS_samples_FINAL.txt"))

} else {
  corrected_calls <- data$gene_calls
}

if(!is.null(data$focal)){
focalGR <- makeGRangesFromDataFrame(data$focal, keep.extra.columns = T)


clustering <- EvoCluster(corrected_calls,
                         eps_pair, 
                         min_pair, 
                         noise_sd, 
                         focalGR, 
                         plot_dir,
                         data$pt_list,
                         show_plot=T)


write_tsv(clustering$all_clusters, file.path(output_dir, "2.dbscan_results.txt"))
write_tsv(clustering$all_focals_alone, file.path(output_dir, "2.focals_not_in_clusters_results.txt"))
} else {
  clustering <- EvoCluster(corrected_calls, eps_pair, min_pair, noise_sd, focalGR =NULL, plot_dir, pt_list, show_plot=F) 
  write_tsv(clustering$all_clusters, file.path(output_dir, "2.dbscan_results.txt"))
}

state_thresholds <- tibble(
  states =  c("hT","T","AT","A","sA","wt","sD","D","HD","H"),
  limits_low = c(4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2, -Inf),
  limits_high = c(Inf, 4.2, 3.8, 3.2, 2.8, 2.2, 1.8, 1.2, 0.8, 0.2)
)


# Process clusters
cluster_results <- assign_labels(clustering$all_clusters, state_thresholds)

evo_clusters <- evoC_scores(cluster_results, data$pos_evol_tab, data$neg_evol_tab)
write_tsv(evo_clusters, file.path(output_dir, "3.cluster_evoscores.txt"))

# Process focal lesions
focal_results <- assign_labels(clustering$all_focals_alone, state_thresholds)

evo_gene <- evoF_scores(focal_results, data$pos_evol_tab, data$neg_evol_tab)

if(!is.null(data$focal)){
write_tsv(evo_gene, file.path(output_dir, "3.focal_evoscores.txt"))
}

message("Processing complete. Results saved to ", output_dir)

final_scores <- CloEvoScore(evo_clusters,evo_gene, output_dir)
# Save final scores
write_tsv(final_scores, file.path(output_dir, "4.Trajectories.txt"))

message("Processing complete. Results saved to ", output_dir)

