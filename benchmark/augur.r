library(Augur)
library(Seurat)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(viridis)

# load data
data("sc_sim")

# predict with highly variable genes
# augur_h = calculate_auc(sc_sim, classifier = 'rf', subsample_size = 20, n_threads = 4, select_var = F)
# augur_h$AUC

# predict with select variance feature selection
augur_v = calculate_auc(sc_sim, classifier = 'rf', subsample_size = 20, n_threads = 16, select_var = T)
augur_v$AUC

# to visualize the results
plot_lollipop(augur_v)
# plot_scatterplot(augur_h, augur_v, top_n = 0)
# preprocess for plotting umap
sc_sim = FindVariableFeatures(sc_sim)
sc_sim = ScaleData(sc_sim)
sc_sim = RunPCA(sc_sim, features = VariableFeatures(object = sc_sim))
sc_sim <- RunUMAP(object  = sc_sim, dims = 1:5) # delete dims
rownames(sc_sim@meta.data) = names(sc_sim@active.ident)
# plot umap
plot_umap(augur_v, sc_sim, reduction = 'umap')

# differntial prioritization
cat("differntial prioritization start")
# load data
sc = readRDS("/home/icb/zihe.zheng/pertpy_benchmark/data/Bhattacherjee2019.rds")
Idents(sc) = sc@meta.data$label
sc_sub1 = subset(sc, idents = c("Maintenance_Cocaine", "withdraw_15d_Cocaine"))
sc_sub2 = subset(sc, idents = c("Maintenance_Cocaine", "withdraw_48h_Cocaine"))

# compare maintainance vs 15d withdraw
# default
bhattacherjee_results_15 = calculate_auc(sc_sub1, n_threads = 16)
# permute
bhattacherjee_results_15_permute = calculate_auc(sc_sub1, n_threads = 16, augur_mode = 'permute', n_subsamples = 100)
cat("permute 1 done")

# compare maintainance vs 48h withdraw
# default
bhattacherjee_results_48 = calculate_auc(sc_sub2, n_threads = 16)
# permute
bhattacherjee_results_48_permute = calculate_auc(sc_sub2, n_threads = 16, augur_mode = 'permute', n_subsamples = 100)
cat("differntial prioritization done")
# visualize
plot_scatterplot(bhattacherjee_results_15, bhattacherjee_results_48, top_n = 0)
pvals = calculate_differential_prioritization(augur1 = bhattacherjee_results_15, 
                                              augur2 = bhattacherjee_results_48, 
                                              permuted1 = bhattacherjee_results_15_permute, 
                                              permuted2 = bhattacherjee_results_48_permute)
plot_differential_prioritization(pvals)









