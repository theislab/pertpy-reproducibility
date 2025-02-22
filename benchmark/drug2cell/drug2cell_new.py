import time

import blitzgsea as blitz
import pertpy as pt
import scanpy as sc

adata = sc.datasets.pbmc3k_processed()

pt_chembl = pt.md.Drug()
pt_chembl.annotate(adata)

start = time.time()

pt_enricher = pt.tl.Enrichment()
pt_enricher.score(adata)

sc.tl.rank_genes_groups(adata, method="wilcoxon", groupby="louvain")

overrepresentation = pt_enricher.hypergeometric(adata)
enrichment = pt_enricher.gsea(adata)

targets = blitz.enrichr.get_library("GO_Molecular_Function_2021")
targets["MHC class II receptor activity (GO:0032395)"]

pt_enricher.score(adata, targets=targets)
enrichment = pt_enricher.gsea(adata, targets=targets)

runtime = time.time() - start
print(f"Runtime: {runtime:.2f} seconds")
