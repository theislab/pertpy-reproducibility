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

# Load dataset
InstallData(ds = "thp1.eccite")
eccite <- LoadData(ds = "thp1.eccite")
eccite <- UpdateSeuratObject(eccite)

# Preprocessing
# Protein
eccite <- NormalizeData(
  object = eccite,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2)

# RNA
DefaultAssay(object = eccite) <- 'RNA'
eccite <- NormalizeData(object = eccite) %>% FindVariableFeatures() %>% ScaleData()

# Create a random matrix with one row per cell and 50 principal components
set.seed(123)
cells <- Cells(eccite)
rand_mat <- matrix(rnorm(length(cells) * 50, 0, 1), nrow = length(cells), ncol = 50)
rownames(rand_mat) <- cells
eccite[["pca"]] <- CreateDimReducObject(
  embeddings = rand_mat,
  key = "PC_",
  assay = DefaultAssay(eccite)
)

# Mitigating confounding effects
eccite<- CalcPerturbSig(
  object = eccite,
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
DefaultAssay(object = eccite) <- 'PRTB'

# identify cells with no detectable perturbation
eccite <- RunMixscape(
  object = eccite,
  assay = "PRTB",
  slot = "scale.data",
  labels = "gene",
  nt.class.name = "NT",
  min.de.genes = 5,
  iter.num = 10,
  de.assay = "RNA",
  verbose = F,
  prtb.type = "KO")

file.create(snakemake@output[[1]])