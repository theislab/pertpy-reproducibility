import pertpy as pt
import scanpy as sc
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import seaborn as sns

import jax.numpy as jnp
from jax import random
from anndata import AnnData

n_obs = int(snakemake.wildcards.n_obs) // 3

def generate_toy_data(n_guides=3, n_cells_per_group=50):
    dats = []
    for i in range(n_guides):
        # Simulate positive and negative population for each guide
        key = random.PRNGKey(i)
        key1, key2 = random.split(key)

        # Negative first
        poisson_data = random.poisson(key1, lam=0.1, shape=n_cells_per_group * i).astype(jnp.float32)

        # Positive
        gaussian_data = random.normal(key2, shape=n_cells_per_group) * 1.0 + 3
        gaussian_data = gaussian_data.clip(0.0, None)

        # Negative second
        poisson_data_ = random.poisson(key1, lam=0.1, shape=(n_cells_per_group * (n_guides - i - 1))).astype(
            jnp.float32
        )

        # The count vector for one guide is the concatenation of the negative and positive populations
        guide_data = jnp.hstack([poisson_data, gaussian_data, poisson_data_])
        dats.append(guide_data)
    guide_counts = np.array(jnp.vstack(dats)).T

    # Combine Poisson and Gaussian data into one dataset
    adata = AnnData(
        guide_counts,
        obs=pd.DataFrame(index=[f"cell{i+1}" for i in range(guide_counts.shape[0])]),
        var=pd.DataFrame(index=[f"guide{i+1}" for i in range(guide_counts.shape[1])]),
    )
    adata.obs["ground_truth"] = ["guide" + str(i + 1) for i in range(n_guides) for _ in range(n_cells_per_group)]
    return adata

adata = generate_toy_data(n_guides=3, n_cells_per_group=n_obs)
ga = pt.pp.GuideAssignment()
ga.assign_mixture_model(adata)

adata.obs.to_csv(snakemake.output[0])