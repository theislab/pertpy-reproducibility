# pertpy-reproducibility

This repository contains all scripts to reproduce the associated publication.

## Use cases

- `norman` contains all analyses related to "Learning and exploring perturbation representations with pertpy".
- `mcfarland` contains all analyses related to "Pertpy streamlines discovery for complex perturbation experiments"
- `zhang` contains all analyses related to "Pertpy enables deciphering effects of perturbations on cellular systems".

## Reproducibility

The `benchmark` folder has individual tool specific subfolders that each contain a Conda `*_environment.yml` file together with a `*_comparison` notebook where the original implementation is compared against the implementation in pertpy with a suitable metric.

## Run time

The `benchmark` folder further contains all scripts that were run using a Snakemake pipeline to evaluate the runtime requirements of the tools.
In particular, the `visualize_benchmark_results.ipynb` notebook aggregates all results from the `benchmarking_table.csv` which is stored in the tables folder to create the figure.
