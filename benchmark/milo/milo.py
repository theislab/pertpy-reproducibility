import time

import pertpy as pt
import scanpy as sc

adata = pt.dt.stephenson_2021_subsampled()
adata = adata[adata.obs["Status"] != "LPS"].copy()

milo = pt.tl.Milo()
mdata = milo.load(adata)

sc.pp.neighbors(mdata["rna"], use_rep="X_scVI", n_neighbors=150)

start = time.time()

milo.make_nhoods(mdata["rna"], prop=0.1)

mdata = milo.count_nhoods(mdata, sample_col="patient_id")

mdata["rna"].obs["Status"] = mdata["rna"].obs["Status"].cat.reorder_categories(["Healthy", "Covid"])
milo.da_nhoods(mdata, design="~Site+Status", model_contrasts="StatusCovid")

runtime = time.time() - start
print(f"Runtime: {runtime:.2f} seconds")
