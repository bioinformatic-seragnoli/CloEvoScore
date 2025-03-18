
# Purpose: Assign Score to clusters and focal regions based on DBscan cluster position.
# === Helper Functions ===

# Determine state for a given CN value
get_state <- function(value, state_thresholds) {
  state_thresholds %>%
    filter(value >= limits_low & value < limits_high) %>%
    pull(states)
}

# Assign  labels to data
assign_labels <- function(pt_clusters, state_thresholds) {
  # Check if the input data has rows
  if (nrow(pt_clusters) == 0) {
    cat("Error: There are no clusters for the patient.\n")
    return(tibble(pt_name = character(), diagnosis_label = character(), relapse_label = character()))
  }
  
  # Initialize variables
  pts <- unique(pt_clusters$pt_name)
  all_results <- list()
  
  # Iterate over each patient
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
      # Assign labels if clusters are valid
      patient_data <- patient_data %>%
        mutate(
          diagnosis_label = map_chr(diagnosis_noise, ~ get_state(.x, state_thresholds)),
          relapse_label = map_chr(relapse_noise, ~ get_state(.x, state_thresholds))
        )
    }
    
    # Append patient results
    all_results[[pt]] <- patient_data
  }
  
  # Combine results and return
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
      pt_name = unique(clusters$pt_name), # Preserve patient names if they exist
      N_clusters = 0,                     # Set the number of clusters to 0
      Nscore = 0,                     # Assign zero negative pressure
      Pscore = 0                   # Assign zero positive pressure
      
    ))
  }
  
  # Compute POS and NEG scores for each cluster
  clusters <- clusters %>% 
    mutate(
      Pscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        # Match relapse state in pos_tab and retrieve the corresponding diagnosis score
        pos_tab %>% filter(POS == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      }),
      Nscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        # Match relapse state in neg_tab and retrieve the corresponding diagnosis score
        neg_tab %>% filter(NEG == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      })
    )
  
  # Aggregate scores at the cluster level
  clusters_scores <- clusters %>% 
    group_by(pt_name, DBscan_cluster) %>% 
    summarise(
      n_chr = n_distinct(chromosome),    # Count unique chromosomes
      chr_list = str_c(unique(chromosome), collapse = ", "), # List of chromosomes
      genes_x_chr = str_c(n_genes, collapse = ", "),         # List of gene counts per chromosome
      Pscore = unique(Pscore)[1],  # Take the first unique Pscore
      Nscore = unique(Nscore)[1],  # Take the first unique Nscore
      .groups = "drop"
    ) %>% 
    mutate(
      Pscore_chr_weight = Pscore * n_chr,  # Weighted positive score by chromosome count
      Nscore_chr_weight = Nscore * n_chr   # Weighted negative score by chromosome count
    )
  
  # Aggregate scores at the patient level
  pt_score <- clusters_scores %>% 
    group_by(pt_name) %>% 
    summarise(
      pt_name = unique(pt_name),                # Ensure unique patient name
      N_clusters = max(DBscan_cluster, na.rm = TRUE), # Count of clusters for the patient
      Nscore = sum(Nscore_chr_weight, na.rm = TRUE),  # Sum of NEG scores across clusters
      Pscore = sum(Pscore_chr_weight, na.rm = TRUE),  # Sum of POS scores across clusters
    )
  
  return(pt_score) # Return patient-level scores
}



# Compute evolutionary pressure scores
evoF_scores <-  function(genes, pos_tab, neg_tab) {
  if (is.null(genes) || nrow(genes) == 0) {
    message("The genes data frame is empty. Returning an empty data frame.")
    return(tibble(pt_name = character(), Pscore = numeric(), Nscore = numeric())) # Return an empty data frame
  }
  genes <- genes %>%
    mutate(
      Pscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        pos_tab %>% filter(POS == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      }),
      Nscore = map2_dbl(diagnosis_label, relapse_label, ~ {
        neg_tab %>% filter(NEG == .y) %>% pull(.x) %>% as.numeric() %>% coalesce(0)
      })
    )
  
  genes$changement <- ifelse(rowSums(genes[,c("Nscore", "Pscore")]) !=0, 1,0)
  
  genes <- genes %>% 
    rowwise() %>%
    mutate(ref_conc = 
             grepl(CNA, diagnosis_label) | 
             grepl(CNA, relapse_label))
  
  
  filtered_genes <- genes %>% filter(changement == 1 & ref_conc == T)
  
  
  final_genes <- filtered_genes %>%
    group_by(pt_name) %>%
    summarise(
      genes = paste(unique(gene_names), collapse = ", "),
      Nscore = sum(Nscore, na.rm = TRUE),
      Pscore = sum(Pscore, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(final_genes)
}
