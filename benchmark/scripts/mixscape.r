# Load packages.
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)

# load dataset
InstallData(ds = "thp1.eccite")
eccite <- LoadData(ds = "thp1.eccite")

# preprocessing
# protein.
eccite <- NormalizeData(
  object = eccite,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2)

# RNA
DefaultAssay(object = eccite) <- 'RNA'
eccite <- NormalizeData(object = eccite) %>% FindVariableFeatures() %>% ScaleData()

# gene expression-based cell clustering umap
eccite <- RunPCA(object = eccite)
eccite <- RunUMAP(object = eccite, dims = 1:40)

# mitigating confounding effects
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
# Normalize data, find variable features and center data.
DefaultAssay(object = eccite) <- 'PRTB'

# Use variable features from RNA assay.
VariableFeatures(object = eccite) <- VariableFeatures(object = eccite[["RNA"]])
eccite <- ScaleData(object = eccite, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
eccite <- RunPCA(object = eccite, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
eccite <- RunUMAP(
  object = eccite,
  dims = 1:40,
  reduction = 'prtbpca',
  reduction.key = 'prtbumap',
  reduction.name = 'prtbumap')

# identiy cells with no detectable perturbation
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

# Remove non-perturbed cells and run LDA to reduce the dimensionality of the data.
Idents(eccite) <- "mixscape_class.global"
sub <- subset(eccite, idents = c("KO", "NT"))

# Run LDA.
sub <- MixscapeLDA(
  object = sub,
  assay = "RNA",
  pc.assay = "PRTB",
  labels = "gene",
  nt.label = "NT",
  npcs = 10,
  logfc.threshold = 0.25,
  verbose = F)
