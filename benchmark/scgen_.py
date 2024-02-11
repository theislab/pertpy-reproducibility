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
print("training completed")
# save the model
# model.save("saved_models/model_perturbation_prediction.pt", overwrite=True)

# obtain latent representation and show in UMAP
latent_X = model.get_latent_representation()
latent_adata = sc.AnnData(X=latent_X, obs=train_new.obs.copy())
sc.pp.neighbors(latent_adata)
sc.tl.umap(latent_adata)
sc.pl.umap(
    latent_adata,
    color=["label", "cell_type"],
    wspace=0.4,
    frameon=False,
    # save="latentspace_batch32_klw000005_z100__100e.pdf",
)
print('umap completed')
# prediction
pred, delta = model.predict(ctrl_key="ctrl", stim_key="stim", celltype_to_predict="CD4 T cells")
pred.obs["label"] = "pred"

## evaluation of prediction
ctrl_adata = train[((train.obs["cell_type"] == "CD4 T cells") & (train.obs["label"] == "ctrl"))]
stim_adata = train[((train.obs["cell_type"] == "CD4 T cells") & (train.obs["label"] == "stim"))]
eval_adata = ctrl_adata.concatenate(stim_adata, pred)

# embedding all real and predicted cells in one PCA plot
sc.tl.pca(eval_adata)
sc.pl.pca(
    eval_adata,
    color="label",
    frameon=False,
    # save="pred_stim_b32_klw000005_z100__100e.pdf",
)

# mean correlation plot
CD4T = train[train.obs["cell_type"] == "CD4 T cells"]
sc.tl.rank_genes_groups(CD4T, groupby="label", method="wilcoxon")
diff_genes = CD4T.uns["rank_genes_groups"]["names"]["stim"]
print(diff_genes)

condition_key = model.adata_manager.get_state_registry(REGISTRY_KEYS.BATCH_KEY).original_key

r2_value = pt.pl.scg.reg_mean_plot(
    eval_adata,
    condition_key=condition_key,
    axis_keys={"x": "pred", "y": "stim"},
    gene_list=diff_genes[:10],
    labels={"x": "predicted", "y": "ground truth"},
    # path_to_save="./reg_mean.pdf",
    show=True,
    legend=False,
)

r2_value, r_value_diff = pt.pl.scg.reg_mean_plot(
    eval_adata,
    condition_key=condition_key,
    axis_keys={"x": "pred", "y": "stim"},
    gene_list=diff_genes[:10],
    top_100_genes=diff_genes,
    x_coeff=0.91,
    y_coeff=0.76,
    labels={"x": "predicted", "y": "ground truth"},
    # path_to_save="./reg_mean_diff_genes.pdf",
    show=True,
    legend=False,
)

# violin plot for a specific gene
sc.pl.violin(eval_adata, keys="ISG15", groupby="label")



