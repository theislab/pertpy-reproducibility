import time
start = time.time()
from pathlib import Path
from pertpy.tools import Mixscape
from scanpy import AnnData, read_h5ad

import cProfile
import io
import pstats
import time
print("Importing libraries took: ", time.time() - start)

# I/O
if "snakemake" in locals():
    input = snakemake.input[0]
    output = snakemake.output[0]
else:
    input = None
    n_obs = None

profiler = cProfile.Profile()
profiler.enable()

start_data_read = time.time()
adata = read_h5ad(input)
print("Time until data was read: ", time.time() - start_data_read)

# Mitigating confounding effects
start_mixscape_init = time.time()
mixscape_identifier = Mixscape()
print("Time until mixscape was initialized: ", time.time() - start_mixscape_init)

start_perturbation_signature = time.time()
mixscape_identifier.perturbation_signature(
    adata, 
    "perturbation", 
    "NT", 
    split_by="replicate", 
    n_neighbors=20, 
    n_dims=40,
)
print("Time taken until perturbation signature done: ", time.time() - start_perturbation_signature)

start_mixscape = time.time()

# Identify cells with no detectable perturbation
mixscape_identifier.mixscape(
    adata=adata, 
    control="NT", 
    labels="gene_target", 
    layer="X_pert"
)
print("Time taken until mixscape done: ", time.time() - start_mixscape)

# Visualizing perturbation responses with Linear Discriminant Analysis (LDA)
start_lda = time.time()
mixscape_identifier.lda(
    adata=adata, control="NT", labels="gene_target"
)
print("Time taken until lda done: ", time.time() - start_lda)
profiler.disable()
s = io.StringIO()
ps = pstats.Stats(profiler, stream=s).sort_stats("cumtime")
ps.dump_stats(filename=str(snakemake.wildcards.n_obs)+"_mixscape.prof")

print("Time taken entire script: ", time.time() - start)

if output:
    Path(output).touch()
