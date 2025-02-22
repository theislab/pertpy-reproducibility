suppressPackageStartupMessages({
  library(Augur)
  library(Seurat)
  library(tidyverse)
  library(magrittr)
  library(ggplot2)
  library(viridis)
})

# I/O
out_sim <- snakemake@output$out_sim
target_n <- as.integer(snakemake@wildcards$n_sample)

# load data
data("sc_sim")
if (target_n == 0) {
  sc_sim_upsampled <- sc_sim
} else {
  # Upsample
  counts <- GetAssayData(sc_sim, slot = "counts")
  orig_meta <- sc_sim[[]]
  rownames(orig_meta) <- colnames(counts)
  cell_names <- colnames(counts)
  
  # Sample cell names with replacement
  new_cell_names <- sample(cell_names, target_n, replace = TRUE)
  
  # Subset counts matrix and metadata using the sampled cell names
  new_counts <- counts[, new_cell_names]
  upsampled_meta <- orig_meta[new_cell_names, ]
  
  # Modify cell names to ensure they are unique for the columns (cells)
  new_names <- paste0(new_cell_names, "_dup", seq_along(new_cell_names))
  colnames(new_counts) <- new_names
  rownames(upsampled_meta) <- new_names
  
  # Create a new Seurat object with the upsampled data and metadata
  sc_sim_upsampled <- CreateSeuratObject(new_counts, meta.data = upsampled_meta)
}

# predict with select variance feature selection
augur = calculate_auc(sc_sim_upsampled, 
                        classifier = 'rf', 
                        subsample_size = 20, 
                        n_threads = snakemake@threads, 
                        select_var = T,
                        augur_mode='default'
                        )
write.csv(augur$AUC, out_sim)
