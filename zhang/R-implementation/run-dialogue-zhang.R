library(Seurat)
library(DIALOGUE)
library(readr)

# import the zhang_pca csv file and set the first column as the row names 
zhang_pca <- read.csv("~/Documents/GitHub/pertpy-reproducibility/zhang/data/zhang_pca.csv")

# set the first column to be the rownames
rownames(zhang_pca) <- zhang_pca$Cell.barcode
# remove column Cell.barcode
zhang_pca <- zhang_pca[,-1]


zhang_X <- read.csv("~/Documents/GitHub/pertpy-reproducibility/zhang/data/zhang_X.csv")
rownames(zhang_X) <- zhang_X$Cell.barcode
# remove column Cell.barcode
zhang_X <- zhang_X[,-1]



zhang_obs <- read.csv("~/Documents/GitHub/pertpy-reproducibility/zhang/data/zhang_obs.csv")
rownames(zhang_obs) <- zhang_obs$Cell.barcode
# remove column Cell.barcode
zhang_obs <- zhang_obs[,-1]


celltypes <- c('t_Bmem-CD27', 't_CD4_Tcm-LMNA', 't_CD4_Treg-FOXP3', 't_CD8_MAIT-KLRB1',
               't_CD8_Tem-GZMK', 't_CD8_Trm-ZNF683', 't_Tn-LEF1', 't_mono-FCN1',
               't_pB-IGHG1')
rA <- list()

for (ct in celltypes) {
  print(ct)
  # extract cell barcodes with Cluster equal to ct 
  mini_obs <- zhang_obs[zhang_obs$Cluster == ct,]
  mini_X <- zhang_X[rownames(mini_obs),]
  mini_pca <- zhang_pca[rownames(mini_obs),]
  mini_celltype <- make.cell.type(name = ct,
                                  as.matrix(t(mini_X)), #tpm expression matrix
                                  as.character(mini_obs$Sample), # samples list of strings?
                                  X = as.matrix(mini_pca),
                                  metadata = mini_obs,
                                  cellQ = as.vector(mini_obs$n_counts))
  rA[[ct]] <- mini_celltype
}
# create the DIALOGUE cell type objects by filtering in a for loop
param <- DLG.get.param(k = 10, 
                       results.dir = "R_dlg_output/", 
                       conf = c("cellQ") # Confounding factors
                       ) # Phenotype (optional)

R<-DIALOGUE.run(rA = rA, # list of cell.type objects
                main = "DLG_run",
                param = param)


# run dialogue with the same cell type ordering used in the original script 
