# === Purpose: Compute evolutionary trajectories and total positive/negative scores for each patient ===
# This function combines scores from clusters and focal regions, then classifies the trajectory of evolution.

# === Helper Function ===

# Function to classify the evolutionary trajectory based on combined scores
calculate_trajectory <- function(all_scores, med_neg = 1, med_pos = 1) {
  all_scores %>%
    mutate(
      Trajectory = case_when(
        # Case 1: Low Pscore and low Nscore → Stable evolution
        Tot_Pscore <= med_pos & Tot_Nscore <= med_neg ~ "stable",
        # Case 2: High Pscore, low or moderate Nscore OR
        # Both high but Pscore > Nscore → Linear-like evolution
        Tot_Pscore > med_pos & Tot_Nscore <= med_neg | 
          Tot_Pscore > med_pos & Tot_Nscore > med_neg & Tot_Pscore > Tot_Nscore ~ "linear-like",
        
        # Case 3: High Nscore and low/moderate Pscore OR
        # Both high but Nscore ≥ Pscore → Drift-like evolution
        Tot_Pscore <= med_pos & Tot_Nscore > med_neg | Tot_Pscore > med_pos & Tot_Nscore > med_neg & Tot_Pscore <= Tot_Nscore ~ "drift-like"
      )
    )
}

# === Main Function to Combine Scores and Compute Trajectories ===
CloEvoScore <- function(evo_clusters,evo_gene,output_dir) {
  # If gene-level score table has data (i.e., contains more than 1 column)
  if(ncol(evo_gene) > 1) {
    
    # Merge cluster and gene scores by patient name
    combined_scores <- full_join(evo_clusters, evo_gene, by = "pt_name", suffix = c("_Clusters", "_Genes"))
  
    # Replace any missing gene scores (NA) with zeroes, then compute total scores
    combined_scores <- combined_scores %>%
    mutate(
      Pscore_Genes = coalesce(Pscore_Genes, 0), 
      Nscore_Genes = coalesce(Nscore_Genes, 0), 
      Tot_Pscore = Pscore_Clusters + Pscore_Genes,
      Tot_Nscore = Nscore_Clusters + Nscore_Genes)
    
    # Classify patient trajectory based on combined scores
    final_scores <- calculate_trajectory(combined_scores)
  } else {
    
    # Case where no gene-level data is available: use only cluster scores
    combined_scores <- evo_clusters %>%  
      mutate (
      Tot_Pscore = Pscore,
      Tot_Nscore = Nscore)
    
    # Classify patient trajectory based only on cluster-level scores
    final_scores <- calculate_trajectory(combined_scores)
    
  }
  
  
  message("Analysis complete")
  return(final_scores) # Return final scores and trajectory per patient
  
}

