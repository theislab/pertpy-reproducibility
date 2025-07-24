library('DIALOGUE')

# I/O
rA <- readRDS(snakemake@input$rA)

# Subsample cells
n_cells <- as.integer(snakemake@wildcards$n_obs)
set.seed(123)
for (i in 1:length(rA)) {
    cell_indices <- sample(seq_along(rA[[i]]@cells), n_cells, replace = TRUE)
    temp = rA[[i]]
    temp@cells <- rA[[i]]@cells[cell_indices]
    temp@cellQ <- rA[[i]]@cellQ[cell_indices]
    temp@tpm <- rA[[i]]@tpm[, cell_indices]
    temp@X <- rA[[i]]@X[cell_indices, ]
    temp@samples <- rA[[i]]@samples[cell_indices]
    temp@metadata <- rA[[i]]@metadata[cell_indices, ]
    rA[[i]] <- temp
}

# Create the directory if it doesn't exist
path <- paste0("/scratch/peidli/pertpy_benchmark/dialogue/ncells_", snakemake@wildcards$n_obs)
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
}
param <- DLG.get.param(k = 3, conf = c("cellQ"), results.dir = path)
res <- DIALOGUE.run(rA = rA, main = "toyExample", param = param)

file.create(snakemake@output[[1]])