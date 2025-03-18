# Purpose: Correct Copy Number Variations (CNVs) and prepare the final dataset

# === Functions ===

correct_cn <- function(gene_calls, pp) {
  # Convert to data.table for better performance
  setDT(gene_calls)
  setDT(pp)
  
  # Get valid sample names
  samples_valid <- setdiff(names(gene_calls), names(gene_calls)[1:4])
  
  # Preallocazione della memoria per il risultato
  corrected_calls <- copy(gene_calls)
  
  # Estrarre pazienti e fasi in un'unica operazione
  patient_info <- data.table(
    sample = samples_valid,
    phase = str_extract(samples_valid, "[EeRrDd]"),
    patient = str_remove(samples_valid, "_E|_R|_D")
  )
  
  # Determinare colonne di purezza e ploidia
  patient_info[, purity_col := ifelse(phase %in% c("E", "e", "d", "D"), "Purity_D", "Purity_R")]
  patient_info[, diploid_col := ifelse(phase %in% c("E", "e", "d", "D"), "Ploidy_D", "Ploidy_R")]
  
  # Unire i dati con pp per ottenere purezza e ploidia
  patient_info <- merge(patient_info, pp, by.x = "patient", by.y = "pt_name", all.x = TRUE)
  
  # Normalizzare i valori
  for (row in 1:nrow(patient_info)) {
    sample_name <- patient_info$sample[row]
    purity <- patient_info[row, get(purity_col)] / 100
    ploidy <- patient_info[row, get(diploid_col)] / 100
    
    # Calcolare direttamente sulla colonna senza creare nuovi oggetti
    corrected_calls[[sample_name]] <- (((gene_calls[[sample_name]] - 2) / purity) + 2) + ploidy
  }
  
  # Impostare i valori negativi a zero
  corrected_calls[, (samples_valid) := lapply(.SD, function(x) fifelse(x < 0, 0, x)), .SDcols = samples_valid]
  
  # Arrotondare i valori
  corrected_calls[, (samples_valid) := lapply(.SD, round, 3), .SDcols = samples_valid]
  
  message("Correction complete")
  
  return(corrected_calls)
}


