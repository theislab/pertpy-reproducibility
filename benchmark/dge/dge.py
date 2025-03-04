import time
import warnings

import decoupler as dc
import pertpy as pt
import scanpy as sc

warnings.filterwarnings("ignore")

adata = pt.dt.zhang_2021()
# sc.pp.sample(adata, n=100000)

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
