# Purpose: Correct Copy Number Variations (CNVs) and prepare the final dataset

# === Functions ===

# Correct Copy Number Variations
correct_cn <- function(gene_calls, pp) {
  # Get valid sample names
  samples_valid <- names(gene_calls)[-c(1:4)] %>%
    as.vector() %>%
    sort()
  
  # Apply corrections
  corrected_calls <- gene_calls
  
  for (sample in samples_valid) {
    phase <- str_extract(sample, "[EeRrDd]")
    patient <- str_remove(sample, "_E|_R")
    
    purity_col <- ifelse(phase %in% c("[Ee]", "[Dd]"), "Purity_D", "Purity_R")
    diploid_col <- ifelse(phase %in% c("[Ee]", "[Dd]"), "Ploidy_D", "Ploidy_R")
    
    purity <- pp %>%
      filter(pt_name == patient) %>%
      pull(!!sym(purity_col)) / 100
    ploidy <- pp %>%
      filter(pt_name == patient) %>%
      pull(!!sym(diploid_col)) / 100
    
    corrected_calls[[sample]] <- (((gene_calls[[sample]] - 2) / purity) + 2) + ploidy
  }
  
  corrected_calls[corrected_calls < 0] <- 0
  
  rounded_data <- corrected_calls[,-c(1:4)] %>%
    mutate(across(everything(), ~ round(.x, 3)))
  
  corrected_calls <- bind_cols(corrected_calls %>% select(1:4), rounded_data)
  message("Correction complete")
  
  return(corrected_calls)
}


