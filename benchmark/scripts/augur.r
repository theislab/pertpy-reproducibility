library(Augur)
library(Seurat)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(viridis)

# I/O
Bhattacherjee2019.rds <- "/home/icb/zihe.zheng/pertpy_benchmark/data/Bhattacherjee2019.rds"
out_cocaine <- "../results/augur_cocaine.csv"
out_sim <- "../results/augur_sim.csv"

# load data
data("sc_sim")

# predict with select variance feature selection
augur_v = calculate_auc(sc_sim, classifier = 'rf', subsample_size = 20, n_threads = 16, select_var = T)
write.csv(augur_v$AUC, out_sim)

# differntial prioritization
# load data
sc = readRDS(Bhattacherjee2019.rds)
Idents(sc) = sc@meta.data$label
sc_sub1 = subset(sc, idents = c("Maintenance_Cocaine", "withdraw_15d_Cocaine"))
sc_sub2 = subset(sc, idents = c("Maintenance_Cocaine", "withdraw_48h_Cocaine"))

# compare maintainance vs 15d withdraw
# default
bhattacherjee_results_15 = calculate_auc(sc_sub1, n_threads = 16)
# permute
bhattacherjee_results_15_permute = calculate_auc(sc_sub1, n_threads = 16, augur_mode = 'permute', n_subsamples = 100)

# compare maintainance vs 48h withdraw
# default
bhattacherjee_results_48 = calculate_auc(sc_sub2, n_threads = 16)
# permute
bhattacherjee_results_48_permute = calculate_auc(sc_sub2, n_threads = 16, augur_mode = 'permute', n_subsamples = 100)

pvals = calculate_differential_prioritization(augur1 = bhattacherjee_results_15,
                                              augur2 = bhattacherjee_results_48,
                                              permuted1 = bhattacherjee_results_15_permute,
                                              permuted2 = bhattacherjee_results_48_permute)

write.csv(pvals, out_cocaine)