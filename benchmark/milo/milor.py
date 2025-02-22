import time

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pertpy as pt
import rpy2.robjects as robjects
import scanpy as sc
import seaborn as sns
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from scipy import io, sparse

# Initialize R packages
base = importr("base")
SingleCellExperiment = importr("SingleCellExperiment")
miloR = importr("miloR")
dplyr = importr("dplyr")
ggplot2 = importr("ggplot2")

pandas2ri.activate()

adata = pt.dt.stephenson_2021_subsampled()
adata = adata[adata.obs["Status"] != "LPS"].copy()

# Convert to R objects
obs = pandas2ri.py2rpy(adata.obs)
obsm_scvi = pandas2ri.py2rpy(adata.obsm["X_scVI"])
obsm_umap = pandas2ri.py2rpy(adata.obsm["X_umap"])

# R code through rpy2
r_code = """
sce <- SingleCellExperiment(
    colData = obs,
    reducedDims = list(
        scvi = obsm_scvi,
        umap = obsm_umap
    )
)
saveRDS(sce, 'data/stephenson_2021_subsampled.rds')

sce <- readRDS('data/stephenson_2021_subsampled.rds')
counts(sce) <- NULL
logcounts(sce) <- NULL
milo.obj <- Milo(sce)

start_time <- Sys.time()

milo.obj <- buildGraph(milo.obj, k=150, d=10, reduced.dim="scvi")
milo.obj <- makeNhoods(milo.obj, k=150, d=10, refined=TRUE, prop=0.1, refinement_scheme="graph")

patient_metadata <- data.frame(colData(sce)) %>%
    dplyr::select(patient_id, Status, Site)
milo.obj <- countCells(milo.obj, samples="patient_id", meta.data=patient_metadata)

design_df <- patient_metadata %>% distinct()
rownames(design_df) <- design_df$patient_id
design_df$Status <- factor(design_df$Status, levels=c('Healthy', 'Covid'))

milo_res <- testNhoods(
    milo.obj,
    design=~Site+Status,
    design.df=design_df,
    reduced.dim="scvi",
    fdr.weighting="graph-overlap"
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units='mins')
print(paste('Total time elapsed:', round(elapsed, 2), 'minutes'))
"""

robjects.r(r_code)
