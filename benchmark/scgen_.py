### run in env scgen

import pertpy as pt
import scanpy as sc
from scvi import REGISTRY_KEYS

# load data
train = pt.dt.kang_2018()
train.raw = train.copy()
sc.pp.normalize_total(train)

# remove simulated CD4 cells for tutorial
train_new = train[~((train.obs["cell_type"] == "CD4 T cells") & (train.obs["label"] == "stim"))]
train_new = train_new.copy()

# preprocessing data
pt.tl.SCGEN.setup_anndata(train_new, batch_key="label", labels_key="cell_type")

# creating the model
model = pt.tl.SCGEN(train_new)

# training the model
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
    accelerator="cpu",
)

# obtain latent representation and show in UMAP
latent_X = model.get_latent_representation()
latent_adata = sc.AnnData(X=latent_X, obs=train_new.obs.copy())
sc.pp.neighbors(latent_adata)
sc.tl.umap(latent_adata)

# prediction
pred, delta = model.predict(ctrl_key="ctrl", stim_key="stim", celltype_to_predict="CD4 T cells")
pred.obs["label"] = "pred"

## evaluation of prediction
ctrl_adata = train[((train.obs["cell_type"] == "CD4 T cells") & (train.obs["label"] == "ctrl"))]
stim_adata = train[((train.obs["cell_type"] == "CD4 T cells") & (train.obs["label"] == "stim"))]
eval_adata = ctrl_adata.concatenate(stim_adata, pred)

# embedding all real and predicted cells in one PCA plot
sc.tl.pca(eval_adata)

# mean correlation plot
CD4T = train[train.obs["cell_type"] == "CD4 T cells"]
sc.tl.rank_genes_groups(CD4T, groupby="label", method="wilcoxon")
diff_genes = CD4T.uns["rank_genes_groups"]["names"]["stim"]
print(diff_genes)