import muon as mu
import pertpy as pt
import scanpy as sc

# Load dataset
mdata = pt.dt.papalexi_2021()

# Preprocessing
# RNA
sc.pp.highly_variable_genes(mdata["rna"], n_top_genes=2000, flavor='seurat_v3', subset=True)
sc.pp.normalize_total(mdata["rna"])
sc.pp.log1p(mdata["rna"])
mdata["rna"].layers["scaled"] = mdata["rna"].X.copy()
sc.pp.scale(mdata["rna"], layer="scaled")

# Protein
mu.prot.pp.clr(mdata["adt"])

# Gene expression-based cell clustering UMAP
sc.pp.pca(mdata["rna"], n_comps=50, layer="scaled")
sc.pp.neighbors(mdata["rna"], metric="cosine")
sc.tl.umap(mdata["rna"])

# Mitigating confounding effects
mixscape_identifier = pt.tl.Mixscape()
mixscape_identifier.perturbation_signature(
    mdata["rna"], "perturbation", "NT", split_by="replicate", n_neighbors=20, n_dims=40,
)

adata_pert = mdata["rna"].copy()
adata_pert.X = adata_pert.layers["X_pert"]
sc.pp.pca(adata_pert)
sc.pp.neighbors(adata_pert, metric="cosine")
sc.tl.umap(adata_pert)

# Identify cells with no detectable perturbation
mixscape_identifier.mixscape(
    adata=mdata["rna"], control="NT", labels="gene_target", layer="X_pert", split_by="replicate"
)

# Visualizing perturbation responses with Linear Discriminant Analysis (LDA)
mixscape_identifier.lda(
    adata=mdata["rna"], control="NT", labels="gene_target", split_by="replicate"
)
