from pathlib import Path
from pertpy.tl import Mixscape
from scanpy import AnnData, read_h5ad

# I/O
if "snakemake" in locals():
    input = snakemake.input[0]
    output = snakemake.output[0]
else:
    input = None
    n_obs = None

adata = read_h5ad(input)

# Mitigating confounding effects
mixscape_identifier = Mixscape()
mixscape_identifier.perturbation_signature(
    adata, "perturbation", "NT", split_by="replicate", n_neighbors=20, n_dims=40,
)

adata_pert = adata.copy()
adata_pert.X = adata_pert.layers["X_pert"]

# Identify cells with no detectable perturbation
mixscape_identifier.mixscape(
    adata=adata, control="NT", labels="gene_target", layer="X_pert"
)

if output:
    Path(output).touch()