import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import mudata as mu
import pertpy as pt
import arviz as az
import pandas as pd

# load data
haber_cells = pt.dt.haber_2017_regions()

# set up model
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(
    haber_cells,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_label",
    sample_identifier="batch",
    covariate_obs=["condition"],
)

# Select control and salmonella data
sccoda_data.mod["coda_salm"] = sccoda_data["coda"][
    sccoda_data["coda"].obs["condition"].isin(["Control", "Salmonella"])
].copy()

sccoda_data = sccoda_model.prepare(
    sccoda_data,
    modality_key="coda_salm",
    formula="condition",
    reference_cell_type="Goblet",
)

# Run MCMC
sccoda_model.run_nuts(sccoda_data, modality_key="coda_salm")

# show results
sccoda_model.summary(sccoda_data, modality_key="coda_salm")
sccoda_model.credible_effects(sccoda_data, modality_key="coda_salm")
# pt.pl.coda.effects_barplot(sccoda_data, modality_key="coda_salm", parameter="Final Parameter")
sccoda_model.set_fdr(sccoda_data, modality_key="coda_salm", est_fdr=0.4)
sccoda_model.summary(sccoda_data, modality_key="coda_salm")

# saving result
# path = "test"
# sccoda_data.write_h5mu(path)
# loading
# sccoda_data_2 = mu.read_h5mu(path)
# sccoda_model.summary(sccoda_data_2, modality_key="coda_salm")

# test_model = pt.tl.Sccoda()
# test_model.get_intercept_df(sccoda_data_2, modality_key="coda_salm")
# test_model.get_effect_df(sccoda_data, modality_key="coda_salm")
# sccoda_data["coda_salm"].varm["intercept_df"]

######################################### SECOND SCRIPT ##########################################

# Load data
haber_cells = pt.dt.haber_2017_regions()

# Convert data to mudata object
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(
    haber_cells,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_label",
    sample_identifier="batch",
    covariate_obs=["condition"],
)

# Select control and salmonella data as one modality
sccoda_data.mod["coda_salm"] = sccoda_data["coda"][
    sccoda_data["coda"].obs["condition"].isin(["Control", "Salmonella"])
].copy()
# print(sccoda_data)

# plot boxplots
# pt.pl.coda.boxplots(sccoda_data, modality_key="coda", feature_name="condition", add_dots=True)
# plt.show()

# Stacked barplot for each sample
# pt.pl.coda.stacked_barplot(sccoda_data, modality_key="coda", feature_name="samples")
# plt.show()

# Stacked barplot for the levels of "Condition"
# pt.pl.coda.stacked_barplot(sccoda_data, modality_key="coda", feature_name="condition")
# plt.show()

# model all three diseases at once
sccoda_data = sccoda_model.prepare(
    sccoda_data,
    modality_key="coda",
    formula="condition",
    reference_cell_type="Endocrine",
)
sccoda_model.run_nuts(sccoda_data, modality_key="coda")
sccoda_model.summary(sccoda_data, modality_key="coda")

# Set salmonella infection as "default" category
sccoda_data = sccoda_model.prepare(
    sccoda_data,
    modality_key="coda_salm",
    formula="C(condition, Treatment('Salmonella'))",
    reference_cell_type="Goblet",
)
sccoda_model.run_nuts(sccoda_data, modality_key="coda_salm")
sccoda_model.summary(sccoda_data, modality_key="coda_salm")

# switching reference cell type
sccoda_data = sccoda_model.prepare(
    sccoda_data,
    modality_key="coda_salm",
    formula="condition",
    reference_cell_type="Enterocyte",
)
sccoda_model.run_nuts(sccoda_data, modality_key="coda_salm")
sccoda_model.summary(sccoda_data, modality_key="coda_salm", est_fdr=0.4)

## result analysis
sccoda_data = sccoda_model.prepare(
    sccoda_data,
    modality_key="coda_salm",
    formula="condition",
    reference_cell_type="Goblet",
)
sccoda_model.run_nuts(sccoda_data, modality_key="coda_salm")

# extended summary
sccoda_model.summary(sccoda_data, modality_key="coda_salm", hdi_prob=0.8, extended=True)

# diagnostic plotting
# salm_arviz = sccoda_model.make_arviz(sccoda_data, modality_key="coda_salm")
# az.plot_trace(
#     salm_arviz,
#     divergences=False,
#     var_names=["alpha", "beta"],
#     coords={"cell_type": salm_arviz.posterior.coords["cell_type_nb"]},
# )
# plt.tight_layout()
# plt.show()

# Run scCODA with each cell type as the reference
cell_types = sccoda_data["coda_salm"].var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["times_credible"]).fillna(0)

for ct in cell_types:
    print(f"Reference: {ct}")

    # Run inference
    sccoda_data = sccoda_model.prepare(
        sccoda_data,
        modality_key="coda_salm",
        formula="condition",
        reference_cell_type=ct,
    )
    sccoda_model.run_nuts(sccoda_data, modality_key="coda_salm")

    # Select credible effects
    cred_eff = sccoda_model.credible_effects(sccoda_data, modality_key="coda_salm")
    cred_eff.index = cred_eff.index.droplevel(level=0)

    # add up credible effects
    results_cycle["times_credible"] += cred_eff.astype("int")

# Calculate percentages
results_cycle["pct_credible"] = results_cycle["times_credible"] / len(cell_types)
results_cycle["is_credible"] = results_cycle["pct_credible"] > 0.5
print(results_cycle)










