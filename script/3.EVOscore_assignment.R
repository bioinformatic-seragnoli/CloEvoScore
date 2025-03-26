
# === Purpose: Assign scores to DBSCAN-derived clusters and focal regions based on their copy number states ===
# This is used to infer evolutionary pressure (positive/negative) from CNV data.

# === Helper Functions ===

# Determine CNV state (label) based on thresholds
get_state <- function(value, state_thresholds) {
  state_thresholds %>%
    filter(value >= limits_low & value < limits_high) %>%
    pull(states)
}

# Assign diagnosis and relapse CN state labels to clusters
assign_labels <- function(pt_clusters, state_thresholds) {
  # If no clusters are provided
  if (nrow(pt_clusters) == 0) {
    cat("Error: There are no clusters for the patient.\n")
    return(tibble(pt_name = character(), diagnosis_label = character(), relapse_label = character()))
  }
  
  # List of patients to process
  pts <- unique(pt_clusters$pt_name)
  all_results <- list()
  
  # Loop over each patient
  for (pt in pts) {
    message("Processing patient: ", pt)
    patient_data <- pt_clusters %>% filter(pt_name == pt)
    
    # Handle cases where there are no valid clusters for the patient
    if (nrow(patient_data) == 0) {
      patient_data <- tibble(
        pt_name = pt,
        diagnosis_label = "wt",
        relapse_label = "wt"
      )
    } else {
      # Assign state labels based on diagnosis and relapse CN values
      patient_data <- patient_data %>%
        mutate(
          diagnosis_label = map_chr(diagnosis_noise, ~ get_state(.x, state_thresholds)),
          relapse_label = map_chr(relapse_noise, ~ get_state(.x, state_thresholds))
        )
    }
    
    # Store results per patient
    all_results[[pt]] <- patient_data
  }
  
  # Compute evolutionary pressure scores for clusters
  bind_rows(all_results)
}


# Compute evolutionary pressure scores
evoC_scores <- function(clusters, pos_tab, neg_tab) {
  # Filter clusters to exclude noise (DBscan_cluster == 0) and ensure minimum genes threshold
  clusters <- clusters %>% 
    filter(DBscan_cluster != 0, n_genes >= 30)
  
  # Handle the case where no valid clusters remain after filtering
  if (nrow(clusters) == 0) {
    message("No valid clusters found. Assigning zero positive and negative pressure scores.")
    return(tibble(
      pt_name = unique(clusters$pt_name),
      N_clusters = 0,                    
      Nscore = 0,                  
      Pscore = 0                   
    
    ))
  }
  
  # Assign Pscore (positive) and Nscore (negative) per cluster based on diagnosis/relapse state transitions
  clusters <- clusters %>% 
    mutate(
      Pscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        pos_tab %>% filter(POS == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      }),
      Nscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        neg_tab %>% filter(NEG == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      })
    )
  
  # Aggregate scores at the cluster level
  clusters_scores <- clusters %>% 
    group_by(pt_name, DBscan_cluster) %>% 
    summarise(
      n_chr = n_distinct(chromosome),    
      chr_list = str_c(unique(chromosome), collapse = ", "), 
      genes_x_chr = str_c(n_genes, collapse = ", "),         
      Pscore = unique(Pscore)[1],  
      Nscore = unique(Nscore)[1],  
      .groups = "drop"
    ) %>% 
    #  # Weight Pscore and Nscore by number of chromosomes
    mutate(
      Pscore_chr_weight = Pscore * n_chr, 
      Nscore_chr_weight = Nscore * n_chr  
    )
  
  # Aggregate scores at the patient level
  pt_score <- clusters_scores %>% 
    group_by(pt_name) %>% 
    summarise(
      pt_name = unique(pt_name),                
      N_clusters = max(DBscan_cluster, na.rm = TRUE),
      Nscore = sum(Nscore_chr_weight, na.rm = TRUE),
      Pscore = sum(Pscore_chr_weight, na.rm = TRUE),
    )
  
  return(pt_score) # Return patient-level scores
}



# Compute evolutionary pressure scores for unclustered Diseae-related genes
evoF_scores <-  function(genes, pos_tab, neg_tab) {
  #Handle empty input
  if (is.null(genes) || nrow(genes) == 0) {
    message("The genes data frame is empty. Returning an empty data frame.")
    return(tibble(pt_name = character(), Pscore = numeric(), Nscore = numeric())) 
  }
  # Assign scores to each gene based on its state transition
  genes <- genes %>%
    mutate(
      Pscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        pos_tab %>% filter(POS == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      }),
      Nscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        neg_tab %>% filter(NEG == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      })
    )
  # Identify genes that actually underwent a state change
  genes$changement <- ifelse(rowSums(genes[,c("Nscore", "Pscore")]) !=0, 1,0)
  
  genes <- genes %>% 
    rowwise() %>%
    mutate(ref_conc = 
             grepl(CNA, diagnosis_label) | 
             grepl(CNA, relapse_label))
  
  # Filter genes with state changes and consistent CNA pattern
  filtered_genes <- genes %>% filter(changement == 1 & ref_conc == T)
  
  # Summarize results at patient level
  final_genes <- filtered_genes %>%
    group_by(pt_name) %>%
    summarise(
      genes = paste(unique(gene_names), collapse = ", "),
      Nscore = sum(Nscore, na.rm = TRUE),
      Pscore = sum(Pscore, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(final_genes) # Return per-patient focal gene scores
}
