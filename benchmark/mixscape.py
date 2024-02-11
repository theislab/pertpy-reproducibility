import muon as mu
import pertpy as pt
import scanpy as sc

# load dataset
mdata = pt.dt.papalexi_2021()
mdata

# preprocessing
# RNA
mdata["rna"]
sc.pp.normalize_total(mdata["rna"])
sc.pp.log1p(mdata["rna"])
sc.pp.highly_variable_genes(mdata["rna"], subset=True)

# Protein
mdata["adt"]
mu.prot.pp.clr(mdata["adt"])

# gene expression-based cell clustering umap
sc.pp.pca(mdata["rna"])
sc.pp.neighbors(mdata["rna"], metric="cosine")
sc.tl.umap(mdata["rna"])
sc.pl.umap(mdata["rna"], color=["replicate", "Phase", "perturbation"])

# regress out confounding effects
mixscape_identifier = pt.tl.Mixscape()
mixscape_identifier.perturbation_signature(mdata["rna"], "perturbation", "NT", "replicate")
adata_pert = mdata["rna"].copy()
adata_pert.X = adata_pert.layers["X_pert"]
sc.pp.pca(adata_pert)
sc.pp.neighbors(adata_pert, metric="cosine")
sc.tl.umap(adata_pert)
sc.pl.umap(adata_pert, color=["replicate", "Phase", "perturbation"])

# identiy cells with no detectable perturbation
mixscape_identifier.mixscape(adata=mdata["rna"], control="NT", labels="gene_target", layer="X_pert")
pt.pl.ms.barplot(mdata["rna"], guide_rna_column="NT")

# inspecting mixscape results
pt.pl.ms.perturbscore(adata=mdata["rna"], labels="gene_target", target_gene="IFNGR2", color="orange")
pt.pl.ms.violin(
    adata=mdata["rna"],
    target_gene_idents=["NT", "IFNGR2 NP", "IFNGR2 KO"],
    groupby="mixscape_class",
)
pt.pl.ms.heatmap(
    adata=mdata["rna"],
    labels="gene_target",
    target_gene="IFNGR2",
    layer="X_pert",
    control="NT",
)
mdata["adt"].obs["mixscape_class_global"] = mdata["rna"].obs["mixscape_class_global"]
pt.pl.ms.violin(
    adata=mdata["adt"],
    target_gene_idents=["NT", "JAK2", "STAT1", "IFNGR1", "IFNGR2", "IRF1"],
    keys="PDL1",
    groupby="gene_target",
    hue="mixscape_class_global",
)

# visualizing perturbation responses with Linear Discriminant Analysis (LDA)
mixscape_identifier = pt.tl.Mixscape()
mixscape_identifier.lda(adata=mdata["rna"], control="NT", labels="gene_target", layer="X_pert")
pt.pl.ms.lda(adata=mdata["rna"], control="NT")











