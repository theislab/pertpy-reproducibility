import warnings
warnings.filterwarnings("ignore")
import pertpy as pt
import scanpy as sc
import pandas as pd
import numpy as np

# I/O
n_sample = int(snakemake.wildcards.n_sample)
output = snakemake.output[0]

# load data
adata = pt.dt.sc_sim_augur()
if n_sample != 0:
    sc.pp.sample(adata, n=n_sample, rng=0, replace=True)
adata.obs_names_make_unique()
ag_rfc = pt.tl.Augur("random_forest_classifier")
loaded_data = ag_rfc.load(adata)

v_adata, v_results = ag_rfc.predict(loaded_data, 
                                    subsample_size=20, 
                                    n_threads=snakemake.threads, 
                                    augur_mode='default', 
                                    select_variance_features=True)
print(v_results["summary_metrics"])
v_results["summary_metrics"].to_csv(output)