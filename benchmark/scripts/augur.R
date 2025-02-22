library(Augur)
library(Seurat)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(viridis)

# I/O
out_sim <- snakemake@output$out_sim

# load data
data("sc_sim")

# predict with select variance feature selection
augur_v = calculate_auc(sc_sim, classifier = 'rf', subsample_size = 20, 
                        n_threads = snakemake@threads, select_var = T)
write.csv(augur_v$AUC, out_sim)
