import numpy as np
import pandas as pd
import pertpy as pt
import scanpy as sc
from pathlib import Path

# I/O
if "snakemake" in locals():
    output = snakemake.output[0]
    n_obs = int(snakemake.wildcards.n_obs)
else:
    output = None
    n_obs = None

adata = pt.dt.dialogue_example()
adata.obs["dlg_sample"] = ["_".join(x.split(".")[:-1]) for x in adata.obs.index]
adata.obs["dlg_sample"] = adata.obs["dlg_sample"].astype("category")

if n_obs:
    if n_obs < 1e6:
        sc.pp.sample(adata, n=n_obs, rng=0, replace=True)
    else:
        # sample function fails for large n_obs
        idx = np.random.choice(adata.obs.index, n_obs, replace=True)
        adata = adata[idx, :]
adata.obs_names_make_unique()
# Random PCA
adata.obsm['X_pca'] = np.random.normal(size=(adata.n_obs, 30))

# ensure that every cell type is represented in every sample, filter out samples which are missing one
isecs = pd.crosstab(adata.obs["subset"], adata.obs["dlg_sample"])
adata = adata[adata.obs['dlg_sample'].isin(list(isecs.columns[isecs.all(axis=0)]))]

dl = pt.tl.Dialogue(
    sample_id="dlg_sample",
    celltype_key="subset",
    n_counts_key="nCount_RNA",
    n_mpcs=3,
)
adata, mcps, ws, ct_subs = dl.calculate_multifactor_PMD(
    adata, ct_order=["A", "B", "C"], normalize=True
)

if output:
    Path(output).touch()