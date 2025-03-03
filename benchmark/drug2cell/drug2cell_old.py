import time
from pathlib import Path

import blitzgsea as blitz
import drug2cell as d2c
import scanpy as sc

# I/O
if "snakemake" in locals():
    output = snakemake.output[0]
    n_obs = int(snakemake.wildcards.n_obs)
else:
    output = None
    n_obs = None

adata = sc.datasets.pbmc3k_processed()
sc.pp.sample(adata, n=n_obs, rng=0, replace=True)

start = time.time()

d2c.score(adata, use_raw=True)

sc.tl.rank_genes_groups(adata.uns["drug2cell"], method="wilcoxon", groupby="louvain")

targets = blitz.enrichr.get_library("GO_Molecular_Function_2021")
targets["MHC class II receptor activity (GO:0032395)"]

d2c.score(adata, targets=targets, use_raw=True)

sc.tl.rank_genes_groups(adata, method="wilcoxon", groupby="louvain", use_raw=True)

overrepresentation = d2c.hypergeometric(adata)
overrepresentation["B cells"]

del adata.uns["drug2cell"]

enrichment, plot_gsea_args = d2c.gsea(adata, targets=targets)

runtime = time.time() - start
print(f"Runtime: {runtime:.2f} seconds")

if output:
    Path(output).touch()
