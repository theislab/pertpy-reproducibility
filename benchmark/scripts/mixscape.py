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

# protein
mdata["adt"]
mu.prot.pp.clr(mdata["adt"])

# gene expression-based cell clustering umap
sc.pp.pca(mdata["rna"])
sc.pp.neighbors(mdata["rna"], metric="cosine")
sc.tl.umap(mdata["rna"])

# mitigating confounding effects
mixscape_identifier = pt.tl.Mixscape()
mixscape_identifier.perturbation_signature(
    mdata["rna"], "perturbation", "NT", "replicate"
)
adata_pert = mdata["rna"].copy()
adata_pert.X = adata_pert.layers["X_pert"]
sc.pp.pca(adata_pert)
sc.pp.neighbors(adata_pert, metric="cosine")
sc.tl.umap(adata_pert)

# identiy cells with no detectable perturbation
mixscape_identifier.mixscape(
    adata=mdata["rna"], control="NT", labels="gene_target", layer="X_pert"
)

# visualizing perturbation responses with Linear Discriminant Analysis (LDA)
mixscape_identifier = pt.tl.Mixscape()
mixscape_identifier.lda(
    adata=mdata["rna"], control="NT", labels="gene_target", layer="X_pert"
)
