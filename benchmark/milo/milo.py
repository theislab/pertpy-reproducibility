import time
from pathlib import Path

import pertpy as pt
import scanpy as sc

# I/O
n_obs = int(snakemake.wildcards.n_obs)
output = snakemake.output[0]
adata = sc.read(snakemake.input[0])

milo = pt.tl.Milo()
mdata = milo.load(adata)

sc.pp.neighbors(mdata["rna"], use_rep="X_scVI", n_neighbors=150)

milo.make_nhoods(mdata["rna"], prop=0.1)
mdata = milo.count_nhoods(mdata, sample_col="patient_id")
mdata["rna"].obs["Status"] = mdata["rna"].obs["Status"].cat.reorder_categories(["Healthy", "Covid"])
milo.da_nhoods(mdata, design="~Site+Status", model_contrasts="StatusCovid")

if output:
    Path(output).touch()