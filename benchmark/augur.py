import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import pertpy as pt

# load data
adata = pt.dt.sc_sim_augur()
ag_rfc = pt.tl.Augur("random_forest_classifier")
loaded_data = ag_rfc.load(adata)

# predict with highly variable genes
v_adata, v_results = ag_rfc.predict(loaded_data, subsample_size=20, n_threads=16)
print(v_results["summary_metrics"])

# # predict with select variance feature selection
# v_adata, v_results = ag_rfc.predict(
#     loaded_data, subsample_size=20, select_variance_features=True, n_threads=4
# )
# to visualize the results
lollipop = pt.pl.ag.lollipop(v_results)
# scatter = pt.pl.ag.scatterplot(v_results, h_results)
sc.pp.neighbors(adata, use_rep="X")
sc.tl.umap(adata)
sc.pl.umap(adata=v_adata, color="augur_score")

# differntial prioritization
print('differntial prioritization start')
# Load data and create classifier
bhattacherjee_adata = pt.dt.bhattacherjee()
ag_rfc = pt.tl.Augur("random_forest_classifier")

# compare maintainance vs 15d withdraw
# default
bhattacherjee_15 = ag_rfc.load(bhattacherjee_adata, condition_label="Maintenance_Cocaine", treatment_label="withdraw_15d_Cocaine")
bhattacherjee_adata_15, bhattacherjee_results_15 = ag_rfc.predict(bhattacherjee_15, random_state=None, n_threads=16)
print(bhattacherjee_results_15["summary_metrics"].loc["mean_augur_score"].sort_values(ascending=False))

# permute
bhattacherjee_adata_15_permute, bhattacherjee_results_15_permute = ag_rfc.predict(bhattacherjee_15, augur_mode="permute", n_subsamples=100, random_state=None, n_threads=16)

# compare maintainance vs 48h withdraw
# default
bhattacherjee_48 = ag_rfc.load(bhattacherjee_adata, condition_label="Maintenance_Cocaine", treatment_label="withdraw_48h_Cocaine")
bhattacherjee_adata_48, bhattacherjee_results_48 = ag_rfc.predict(bhattacherjee_48, random_state=None, n_threads=16)
print(bhattacherjee_results_48["summary_metrics"].loc["mean_augur_score"].sort_values(ascending=False))

# permute
bhattacherjee_adata_48_permute, bhattacherjee_results_48_permute = ag_rfc.predict(bhattacherjee_48, augur_mode="permute", n_subsamples=100, random_state=None, n_threads=16)

print('differntial prioritization end')
# visualize
scatter = pt.pl.ag.scatterplot(bhattacherjee_results_15, bhattacherjee_results_48)

pvals = ag_rfc.predict_differential_prioritization(
    augur_results1=bhattacherjee_results_15,
    augur_results2=bhattacherjee_results_48,
    permuted_results1=bhattacherjee_results_15_permute,
    permuted_results2=bhattacherjee_results_48_permute,
)

diff = pt.pl.ag.dp_scatter(pvals)