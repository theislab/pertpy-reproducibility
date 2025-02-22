import scanpy as sc
import cinemaot as co
import pertpy as pt
import numpy as np
import time

adata = pt.dt.cinemaot_example()
sc.pp.sample(adata, n=10000, replace=True)
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
