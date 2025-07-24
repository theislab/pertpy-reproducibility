suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
    library(miloR)
    library(Matrix)
})

flag <- snakemake@output[[1]]
target_n <- as.integer(snakemake@wildcards$n_obs)

# Prepare data
data("sim_trajectory", package = "miloR")
## Extract SingleCellExperiment object
traj_sce <- sim_trajectory[['SCE']]
colData(traj_sce) <- DataFrame(sim_trajectory[["meta"]])
colnames(traj_sce) <- colData(traj_sce)$cell_id
redim <- reducedDim(traj_sce, "PCA")
dimnames(redim) <- list(colnames(traj_sce), paste0("PC", c(1:50)))
reducedDim(traj_sce, "PCA") <- redim 

# Upsample
new_indices <- sample(seq_len(ncol(traj_sce)), size = target_n, replace = TRUE)
traj_sce_upsampled <- traj_sce[, new_indices]

# Run Milo
traj_milo <- Milo(traj_sce_upsampled)
traj_milo <- buildGraph(traj_milo, k=150, d=10, reduced.dim="PCA")
traj_milo <- makeNhoods(traj_milo, k=15, d=10, refined=TRUE, prop=0.01, refinement_scheme="graph")
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="Sample")

# Define Design
traj_design <- data.frame(colData(traj_milo))[,c("Sample", "Condition")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$Sample
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
da_results <- testNhoods(traj_milo, design = ~ Condition, design.df = traj_design)

# Do testing
milo_res <- testNhoods(
    traj_milo,
    design=~Condition,
    design.df=traj_design,
    reduced.dim="PCA",
    fdr.weighting="graph-overlap"
)

file.create(flag)