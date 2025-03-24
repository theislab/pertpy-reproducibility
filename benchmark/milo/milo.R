suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
    library(miloR)
    library(Matrix)
})

# Construct SCE object with empry counts and logcounts
obs <- read.csv(snakemake@input$obs)
obsm_scvi <- read.csv(snakemake@input$obsm_scvi)
obsm_umap <- read.csv(snakemake@input$obsm_umap)
sce <- SingleCellExperiment(
    colData = obs,
    reducedDims = list(
        scvi = as.matrix(obsm_scvi),
        umap = as.matrix(obsm_umap)
    )
)
counts(sce) <- NULL
logcounts(sce) <- NULL

milo.obj <- Milo(sce)
milo.obj <- buildGraph(milo.obj, k=150, d=10, reduced.dim="scvi")
milo.obj <- makeNhoods(milo.obj, k=150, d=10, refined=TRUE, prop=0.1, refinement_scheme="graph")

patient_metadata <- data.frame(colData(sce)) %>%
    dplyr::select(patient_id, Status, Site)
milo.obj <- countCells(milo.obj, samples="patient_id", meta.data=patient_metadata)

design_df <- patient_metadata %>% distinct()
rownames(design_df) <- design_df$patient_id
design_df$Status <- factor(design_df$Status, levels=c('Healthy', 'Covid'))

milo_res <- testNhoods(
    milo.obj,
    design=~Site+Status,
    design.df=design_df,
    reduced.dim="scvi",
    fdr.weighting="graph-overlap"
)

write.csv(milo_res, snakemake@output$milo_res)
Matrix::writeMM(nhoods(milo.obj), snakemake@output$milo_res_graph)
