import time
from pathlib import Path

import pertpy as pt
import scanpy as sc

# I/O
if "snakemake" in locals():
    output = snakemake.output[0]
    n_obs = int(snakemake.wildcards.n_obs)
else:
    output = None
    n_obs = 10000

adata = pt.dt.cinemaot_example()
sc.pp.sample(adata, n=n_obs, replace=True)
adata.X = adata.raw.X.copy()

sc.pp.pca(adata)

cot = pt.tl.Cinemaot()
# warm up cache & jit
de = cot.causaleffect(
    adata,
    pert_key="perturbation",
    control="No stimulation",
    return_matching=True,
    thres=1,
    smoothness=3e-5,
    eps=1e-3,
    solver="Sinkhorn",
    preweight_label="cell_type0528",
)

start = time.time()

de = cot.causaleffect(
    adata,
    pert_key="perturbation",
    control="No stimulation",
    return_matching=True,
    thres=1,
    smoothness=3e-5,
    eps=1e-3,
    solver="Sinkhorn",
    preweight_label="cell_type0528",
)

runtime = time.time() - start
print(f"Runtime: {runtime:.2f} seconds")

if output:
    Path(output).touch()