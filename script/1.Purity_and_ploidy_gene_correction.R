# Purpose: Correct Copy Number Variations (CNVs) and prepare the final dataset

# === Functions ===

correct_cn <- function(gene_calls, pp) {
  # Convert to data.table for better performance
  setDT(gene_calls)
  setDT(pp)
  
  # Get valid sample names
  samples_valid <- setdiff(names(gene_calls), names(gene_calls)[1:4])
  
  corrected_calls <- copy(gene_calls)
  
  #Extract patient name and sample phase
  patient_info <- data.table(
    sample = samples_valid,
    phase = str_extract(samples_valid, "[EeRrDd]"),
    patient = str_remove(samples_valid, "_E|_R|_D")
  )
  
  # Find Purity and Ploidy Column Indexes
  patient_info[, purity_col := ifelse(phase %in% c("E", "e", "d", "D"), "Purity_D", "Purity_R")]
  patient_info[, diploid_col := ifelse(phase %in% c("E", "e", "d", "D"), "Ploidy_D", "Ploidy_R")]
  
  # Merge with pp table
  patient_info <- merge(patient_info, pp, by.x = "patient", by.y = "pt_name", all.x = TRUE)
  
  # Normalise values
  for (row in 1:nrow(patient_info)) {
    sample_name <- patient_info$sample[row]
    purity <- patient_info[row, get(purity_col)] / 100
    ploidy <- patient_info[row, get(diploid_col)] / 100
    
    corrected_calls[[sample_name]] <- (((gene_calls[[sample_name]] - 2) / purity) + 2) + ploidy
  }
  
  # Set negative values to 0
  corrected_calls[, (samples_valid) := lapply(.SD, function(x) fifelse(x < 0, 0, x)), .SDcols = samples_valid]
  
  # Rounding values
  corrected_calls[, (samples_valid) := lapply(.SD, round, 3), .SDcols = samples_valid]
  
  message("Correction complete")
  
  return(corrected_calls)
}


