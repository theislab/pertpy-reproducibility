import time
start = time.time()
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from pertpy.tools import Augur
from scanpy import read_h5ad
from anndata import AnnData
from pathlib import Path
print("Importing libraries took: ", time.time() - start)

# I/O
input = snakemake.input[0]
output = snakemake.output[0]

# load data
adata = read_h5ad(input)

ag_rfc = Augur("random_forest_classifier")
loaded_data = ag_rfc.load(adata)

v_adata, v_results = ag_rfc.predict(loaded_data, 
                                    subsample_size=20, 
                                    n_threads=snakemake.threads, 
                                    augur_mode='default', 
                                    select_variance_features=True)
print(v_results["summary_metrics"])

print("Time taken entire script: ", time.time() - start)

if output:
    Path(output).touch()