# Prepare the data for the Milo benchmark, including subsampling and saving the
# data to be read in both R and python.
import pertpy as pt
import scanpy as sc
import pandas as pd
from scipy.io import mmwrite

# I/O
output_h5ad = snakemake.output.h5ad
output_mtx = snakemake.output.mtx
output_obs = snakemake.output.obs
output_obsm_scvi = snakemake.output.obsm_scvi
output_obsm_umap = snakemake.output.obsm_umap
n_obs = int(snakemake.wildcards.n_obs)
print(f"Subsampling {n_obs} cells")

adata = pt.dt.stephenson_2021_subsampled()
adata = adata[adata.obs["Status"] != "LPS"].copy()

# Subsample
if n_obs != 0:
    sc.pp.sample(adata, n=n_obs, rng=0, replace=True)
adata.obs_names_make_unique()
obs = adata.obs.copy()
obsm_scvi = pd.DataFrame(adata.obsm["X_scVI"].copy(), index=adata.obs_names)
obsm_umap = pd.DataFrame(adata.obsm["X_umap"].copy(), index=adata.obs_names)

# Save
mmwrite(output_mtx, adata.X)
adata.write_h5ad(output_h5ad)
obs.to_csv(output_obs)
obsm_scvi.to_csv(output_obsm_scvi)
obsm_umap.to_csv(output_obsm_umap)