import time
import warnings
from pathlib import Path

import decoupler as dc
import pertpy as pt
import scanpy as sc
import numpy as np

warnings.filterwarnings("ignore")

# I/O
if "snakemake" in locals():
    output = snakemake.output[0]
    n_obs = int(snakemake.wildcards.n_obs)
else:
    output = None
    n_obs = None
print(f"n_obs: {n_obs}")

adata = pt.dt.zhang_2021()
if n_obs:
    if n_obs < 3e6:
        sc.pp.sample(adata, n=n_obs, rng=0, replace=True)
    else:
        # sample function fails for large n_obs
        idx = np.random.choice(adata.obs.index, n_obs, replace=True)
        adata = adata[idx, :]

adata = adata[adata.obs["Origin"] == "t", :].copy()
adata = adata[~adata.obs["Patient"].isin(["P010"])]
adata = adata[~adata.obs["Cluster"].isin(["Mix"])]
adata.layers["counts"] = adata.X.copy()

ps = pt.tl.PseudobulkSpace()
pdata = ps.compute(
    adata, target_col="Patient", groups_col="Cluster", layer_key="counts", mode="sum", min_cells=10, min_counts=1000
)

pdata.layers["counts"] = pdata.X.copy()

sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

dc.swap_layer(pdata, "counts", X_layer_key=None, inplace=True)

start = time.time()

pds2 = pt.tl.PyDESeq2(adata=pdata, design="~Efficacy+Treatment")

pds2.fit()

runtime = time.time() - start
print(f"Runtime: {runtime:.2f} seconds")

if output:
    Path(output).touch()
