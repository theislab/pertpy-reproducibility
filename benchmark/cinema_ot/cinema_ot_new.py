import time
start = time.time()
from pathlib import Path
import pertpy as pt
import numpy as np
import scanpy as sc
t = time.time()
print("Importing libraries took: ", t - start)

# I/O
if "snakemake" in locals():
    output = snakemake.output[0]
    n_obs = int(snakemake.wildcards.n_obs)
else:
    output = None
    n_obs = 10000

adata = pt.dt.cinemaot_example()
sc.pp.sample(adata, n=n_obs, replace=True)
adata.obs_names_make_unique()  # throws error with many cells

adata.X = adata.raw.X.copy()

adata.obsm['X_pca'] = np.random.normal(size=(adata.n_obs, 30))
print("Data loaded and preprocessed took: ", time.time() - t)
t = time.time()

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
    preweight_label=None,
    # preweight_label="cell_type0528",
)

start_2 = time.time()

de = cot.causaleffect(
    adata,
    pert_key="perturbation",
    control="No stimulation",
    return_matching=True,
    thres=1,
    smoothness=3e-5,
    eps=1e-3,
    solver="Sinkhorn",
    preweight_label=None,
    # preweight_label="cell_type0528",
)

print("Time taken entire script: ", time.time() - start)
print("Time taken causaleffect minus warmup: ", time.time() - start_2)


if output:
    Path(output).touch()