import time
from pathlib import Path

import cinemaot as co
import numpy as np
import pertpy as pt
import scanpy as sc

# I/O
if snakemake in globals():
    output = snakemake.output[0]
    n_obs = int(snakemake.wildcards.n_obs)
else:
    output = None
    n_obs = 10000

adata = pt.dt.cinemaot_example()
sc.pp.sample(adata, n=n_obs, replace=True)
adata.X = adata.raw.X.copy()

sc.pp.pca(adata)

start = time.time()

cf_orig, ot_orig, de_orig = co.cinemaot.cinemaot_unweighted(
    adata,
    obs_label="perturbation",
    ref_label="IFNb",
    expr_label="No stimulation",
    mode="parametric",
    thres=1,
    smoothness=3e-5,
    eps=1e-3,
    preweight_label="cell_type0528",
)

adata.obsm["cf_orig"] = cf_orig.copy()
adata.obsm["cf_orig"][adata.obs["perturbation"] == "IFNb", :] = np.matmul(
    ot_orig / np.sum(ot_orig, axis=1)[:, None],
    cf_orig[adata.obs["perturbation"] == "No stimulation", :],
)

runtime = time.time() - start
print(f"Runtime: {runtime:.2f} seconds")

if output:
    Path(output).touch()