# Purpose: Compute evolutionary trajectories and final positive/negative scores

# === Helper Functions ===

# Combine results and calculate trajectory
calculate_trajectory <- function(all_scores, med_neg = 1, med_pos = 1) {
  all_scores %>%
    mutate(
      Trajectory = case_when(
        Tot_Pscore <= med_pos & Tot_Nscore <= med_neg ~ "stable",
        Tot_Pscore > med_pos & Tot_Nscore <= med_neg | Tot_Pscore > med_pos & Tot_Nscore > med_neg & Tot_Pscore > Tot_Nscore ~ "linear-like",
        Tot_Pscore <= med_pos & Tot_Nscore > med_neg | Tot_Pscore > med_pos & Tot_Nscore > med_neg & Tot_Pscore <= Tot_Nscore ~ "drift-like"
      )
    )
}

# === Main Processing ===
CloEvoScore <- function(evo_clusters,evo_gene,output_dir) {
  if(ncol(evo_gene) > 1) {
  combined_scores <- full_join(evo_clusters, evo_gene, by = "pt_name", suffix = c("_Clusters", "_Genes"))
  
  combined_scores <- combined_scores %>%
    mutate(
      Pscore_Genes = coalesce(Pscore_Genes, 0), # Replace NA with 0
      Nscore_Genes = coalesce(Nscore_Genes, 0), # Replace NA with 0
      Tot_Pscore = Pscore_Clusters + Pscore_Genes,
      Tot_Nscore = Nscore_Clusters + Nscore_Genes)

  final_scores <- calculate_trajectory(combined_scores)
  } else {
    
    combined_scores <- evo_clusters %>%  
      mutate (
      Tot_Pscore = Pscore,
      Tot_Nscore = Nscore)
    
    final_scores <- calculate_trajectory(combined_scores)
    
  }
  
  
  message("Analysis complete")
  return(final_scores)
  
}

