from pathlib import Path

import muon as mu
import pertpy as pt
import scanpy as sc
import numpy as np

# I/O
if "snakemake" in locals():
    output = snakemake.output[0]
    n_obs = int(snakemake.wildcards.n_obs)
else:
    output = None
    n_obs = None

# Load dataset
mdata = pt.dt.papalexi_2021()
if n_obs:
    from numpy.random import choice
    idx = choice(np.arange(mdata.n_obs), n_obs, replace=True)
    adata = mdata['rna'][idx].copy()
    adata.obs_names_make_unique()
    del mdata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', subset=True)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Randomly generate PCA data
np.random.seed(123)
random_pca = np.random.normal(loc=0, scale=1, size=(adata.n_obs, 50))
adata.obsm["X_pca"] = random_pca

# export
adata.write_h5ad(output)