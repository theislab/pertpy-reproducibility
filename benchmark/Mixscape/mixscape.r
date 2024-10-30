# Load packages
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)

# Load dataset
InstallData(ds = "thp1.eccite")
eccite <- LoadData(ds = "thp1.eccite")

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

# Gene expression-based cell clustering umap
eccite <- RunPCA(object = eccite)
eccite <- RunUMAP(object = eccite, dims = 1:40)

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

# Use variable features from RNA assay
VariableFeatures(object = eccite) <- VariableFeatures(object = eccite[["RNA"]])
eccite <- ScaleData(object = eccite, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data
eccite <- RunPCA(object = eccite, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2D
eccite <- RunUMAP(
  object = eccite,
  dims = 1:40,
  reduction = 'prtbpca',
  reduction.key = 'prtbumap',
  reduction.name = 'prtbumap')

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

# Remove non-perturbed cells
Idents(eccite) <- "mixscape_class.global"
sub <- subset(eccite, idents = c("KO", "NT"))

# run LDA to reduce the dimensionality of the data
sub <- MixscapeLDA(
  object = sub,
  assay = "RNA",
  pc.assay = "PRTB",
  labels = "gene",
  nt.label = "NT",
  npcs = 10,
  logfc.threshold = 0.25,
  verbose = F)

# Save results
SaveH5Seurat(eccite, filename = "output/mixscape_original.h5Seurat")
Convert("output/mixscape_original.h5Seurat", dest = "h5ad")

output_dir <- "gene_csv_files"
dir.create(output_dir, showWarnings = FALSE)
all_genes <- unique(eccite@meta.data$gene)
for (gene in all_genes) {
  # Get the dataframe for the current gene
  df <- Tool(object = eccite, slot = "RunMixscape")[[gene]]
  file_path <- file.path(output_dir, paste0(gene, ".csv"))
  write.csv(df, file = file_path, row.names = TRUE)
}
zip_filename <- "gene_csv_files.zip"
zip(zipfile = zip_filename, files = list.files(output_dir, full.names = TRUE))
