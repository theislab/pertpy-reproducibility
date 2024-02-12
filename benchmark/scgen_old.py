### run in env scgenold

import logging
import scanpy as sc
import scgen

# load data
train = sc.read("./tests/data/train_kang.h5ad",
                backup_url='https://drive.google.com/uc?id=1r87vhoLLq6PXAYdmyyd89zG90eJOFYLk')
train.raw = train.copy()
sc.pp.normalize_total(train)

# remove simulated CD4 cells for tutorial
train_new = train[~((train.obs["cell_type"] == "CD4T") &
                    (train.obs["condition"] == "stimulated"))]
train_new = train_new.copy()

# preprocessing data
scgen.SCGEN.setup_anndata(train_new, batch_key="condition", labels_key="cell_type")

# creating the model
model = scgen.SCGEN(train_new)

# training the model
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25
)

# obtain latent representation and show in UMAP
latent_X = model.get_latent_representation()
latent_adata = sc.AnnData(X=latent_X, obs=train_new.obs.copy())
sc.pp.neighbors(latent_adata)
sc.tl.umap(latent_adata)

# prediction
pred, delta = model.predict(
    ctrl_key='control',
    stim_key='stimulated',
    celltype_to_predict='CD4T'
)
pred.obs['condition'] = 'pred'

## evaluation of prediction
ctrl_adata = train[((train.obs['cell_type'] == 'CD4T') & (train.obs['condition'] == 'control'))]
stim_adata = train[((train.obs['cell_type'] == 'CD4T') & (train.obs['condition'] == 'stimulated'))]
eval_adata = ctrl_adata.concatenate(stim_adata, pred)

# embedding all real and predicted cells in one PCA plot
sc.tl.pca(eval_adata)

# mean correlation plot
CD4T = train[train.obs["cell_type"] =="CD4T"]
sc.tl.rank_genes_groups(CD4T, groupby="condition", method="wilcoxon")
diff_genes = CD4T.uns["rank_genes_groups"]["names"]["stimulated"]
print(diff_genes)