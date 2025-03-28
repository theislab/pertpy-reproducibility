# Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dplyr)
  library(reshape2)
})

# I/O
flag <- snakemake@output[[1]]
target_n <- as.integer(snakemake@wildcards$n_obs)

# Load dataset
InstallData(ds = "thp1.eccite")
eccite <- LoadData(ds = "thp1.eccite")
eccite <- UpdateSeuratObject(eccite)

# Upsample
rna <- GetAssayData(eccite, slot = "counts", assay = "RNA")
adt <- GetAssayData(eccite, slot = "counts", assay = "ADT")
orig_meta <- eccite[[]]
rownames(orig_meta) <- colnames(rna)
cell_names <- colnames(rna)
# Sample cell names with replacement
new_cell_names <- sample(cell_names, target_n, replace = TRUE)
new_rna <- rna[, new_cell_names]
new_adt <- adt[, new_cell_names]
upsampled_meta <- orig_meta[new_cell_names, ]
new_names <- paste0(new_cell_names, "_dup", seq_along(new_cell_names))
colnames(new_rna) <- new_names
colnames(new_adt) <- new_names
rownames(upsampled_meta) <- new_names
eccite_upsampled <- CreateSeuratObject(new_rna, meta.data = upsampled_meta)
# add upsampled data
eccite_upsampled[["RNA"]] <- CreateAssay5Object(counts = new_rna)
eccite_upsampled[["ADT"]] <- CreateAssay5Object(counts = new_adt)

# Preprocessing
# Protein
eccite_upsampled <- Seurat::NormalizeData(
  object = eccite_upsampled,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2)

# RNA
DefaultAssay(object = eccite_upsampled) <- 'RNA'
eccite_upsampled <- NormalizeData(object = eccite_upsampled) %>% FindVariableFeatures() %>% ScaleData()

# Create a random matrix with one row per cell and 50 principal components
set.seed(123)
cells <- Cells(eccite_upsampled)
rand_mat <- matrix(rnorm(length(cells) * 50, 0, 1), nrow = length(cells), ncol = 50)
rownames(rand_mat) <- cells
eccite_upsampled[["pca"]] <- CreateDimReducObject(
  embeddings = rand_mat,
  key = "PC_",
  assay = DefaultAssay(eccite_upsampled)
)

# Mitigating confounding effects
eccite_upsampled<- CalcPerturbSig(
  object = eccite_upsampled,
  assay = "RNA",
  slot = "data",
  gd.class ="gene",
  nt.cell.class = "NT",
  reduction = "pca",
  ndims = 40,
  num.neighbors = 20,
  split.by = "replicate",
  new.assay.name = "PRTB")
# Prepare PRTB assay for dimensionality reduction:
# Normalize data, find variable features and center data
DefaultAssay(object = eccite_upsampled) <- 'PRTB'

# identify cells with no detectable perturbation
eccite_upsampled <- RunMixscape(
  object = eccite_upsampled,
  assay = "PRTB",
  slot = "scale.data",
  labels = "gene",
  nt.class.name = "NT",
  min.de.genes = 5,
  iter.num = 10,
  de.assay = "RNA",
  verbose = F,
  prtb.type = "KO")

file.create(flag)