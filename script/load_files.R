# Purpose: Annotate segments with gene names

# Load data
load_data <- function(gene_tb, pp_tb =NULL, rel_genes = NULL) {
  # Load segment data and gene database from files
  gene_calls <- fread(gene_tb)
  pos_evol_tab <- fread(file.path(data_dir, "POS_pressures.txt"))
  neg_evol_tab <- fread(file.path(data_dir, "NEG_pressures.txt"))
  
  rel_genes <-if(!is.null(rel_genes)) fread(rel_genes) else NULL
  pp <-  if(!is.null(pp_tb)) fread(pp_tb) else NULL
  
  pt_list <- names(gene_calls)[-c(1:4)] %>% str_remove_all(pattern = paste0(c("_[Rr]$", "_[Er]$", "_[Dd]$"), collapse = "|")) %>% unique()
  
  list(gene_calls = gene_calls,
       pos_evol_tab = pos_evol_tab,
       neg_evol_tab = neg_evol_tab, 
       rel_genes = rel_genes, 
       pp=pp,
       pt_list=pt_list)
}



