# === Purpose: Load required files ===
# This function loads segment CNV data, associated pressure score tables, and optional files
# for posterior probabilities and Disease-related genes.

# === Load Data Function ===
load_data <- function(gene_tb, pp_tb =NULL, rel_genes = NULL) {
  # Load the main CNV segment-to-gene mapping table
  gene_calls <- fread(gene_tb)
  
  # Load evolutionary pressure score tables for positive and negative selection
  pos_evol_tab <- fread(file.path(data_dir, "POS_pressures.txt"))
  neg_evol_tab <- fread(file.path(data_dir, "NEG_pressures.txt"))
  
  # Load Disease-related gene table if provided
  rel_genes <-if(!is.null(rel_genes)) fread(rel_genes) else NULL
  pp <-  if(!is.null(pp_tb)) fread(pp_tb) else NULL
  
  # Extract patient names from the CNV table column names
  # Assumes first 4 columns are metadata (e.g., chromosome, start, end, gene)
  # Then removes suffixes like _D, _R, _E from sample names to get unique patient IDs
  pt_list <- names(gene_calls)[-c(1:4)] %>% str_remove_all(pattern = paste0(c("_[Rr]$", "_[Er]$", "_[Dd]$"), collapse = "|")) %>% unique()
  
  # Return all relevant loaded data in a named list
  list(
    gene_calls = gene_calls,         # CNV segments mapped to genes
    pos_evol_tab = pos_evol_tab,     # Positive pressure scoring table
    neg_evol_tab = neg_evol_tab,     # Negative pressure scoring table
    rel_genes = rel_genes,           # Disease-relevant gene table (optional)
    pp = pp,                         # Posterior probability table (optional)
    pt_list = pt_list                # List of patient IDs
  )
}
}



