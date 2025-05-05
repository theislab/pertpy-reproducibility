import time
start = time.time()
from pathlib import Path
from pertpy.tools import Mixscape
from scanpy import AnnData, read_h5ad
print("Importing libraries took: ", time.time() - start)

# I/O
if "snakemake" in locals():
    input = snakemake.input[0]
    output = snakemake.output[0]
else:
    input = None
    n_obs = None

adata = read_h5ad(input)
print("Time until data was read: ", time.time() - start)

# Mitigating confounding effects
mixscape_identifier = Mixscape()
mixscape_identifier.perturbation_signature(
    adata, 
    "perturbation", 
    "NT", 
    split_by="replicate", 
    n_neighbors=20, 
    n_dims=40,
)

print("Time taken until perturbation signature done: ", time.time() - start)
adata_pert = adata.copy()
adata_pert.X = adata_pert.layers["X_pert"]

# Identify cells with no detectable perturbation
mixscape_identifier.mixscape(
    adata=adata, 
    control="NT", 
    labels="gene_target", 
    layer="X_pert"
)
print("Time taken until mixscape done: ", time.time() - start)

# Visualizing perturbation responses with Linear Discriminant Analysis (LDA)
mixscape_identifier.lda(
    adata=adata, control="NT", labels="gene_target"
)

print("Time taken entire script: ", time.time() - start)


if output:
    Path(output).touch()