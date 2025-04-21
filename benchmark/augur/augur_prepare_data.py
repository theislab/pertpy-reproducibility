import warnings
warnings.filterwarnings("ignore")
from pertpy.data import sc_sim_augur
from scanpy.pp import sample
from anndata import AnnData
import pandas as pd

# I/O
n_obs = int(snakemake.wildcards.n_obs)
output = snakemake.output[0]

# load data
adata = sc_sim_augur()
if n_obs != 0:
    sample(adata, n=n_obs, rng=0, replace=True)
adata.obs_names_make_unique()

# export
adata.write_h5ad(output)