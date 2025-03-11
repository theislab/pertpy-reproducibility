import warnings

warnings.filterwarnings("ignore")
import pandas as pd
import pertpy as pt
from sccoda.util import comp_ana as mod

# Load data
haber_cells = pt.dt.haber_2017_regions()
# Extract condition from mouse name and add it as an extra column to the covariates
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(
    haber_cells,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_label",
    sample_identifier="batch",
    covariate_obs=["condition"],
)
sccoda_data.mod["coda_salm"] = sccoda_data["coda"][
    sccoda_data["coda"].obs["condition"].isin(["Control", "Salmonella"])
].copy()

# Select control and salmonella data
data_salm = sccoda_data.mod["coda_salm"]

# run model
model_salm = mod.CompositionalAnalysis(data_salm, formula="condition", reference_cell_type="Goblet")
# run mcmc
sim_results = model_salm.sample_hmc()

# show results
sim_results.summary()
print(sim_results.credible_effects())
sim_results.set_fdr(est_fdr=0.4)
sim_results.summary()

######################################### SECOND SCRIPT ##########################################

# Load data
haber_cells = pt.dt.haber_2017_regions()
# Extract condition from mouse name and add it as an extra column to the covariates
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(
    haber_cells,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_label",
    sample_identifier="batch",
    covariate_obs=["condition"],
)
sccoda_data.mod["coda_salm"] = sccoda_data["coda"][
    sccoda_data["coda"].obs["condition"].isin(["Control", "Salmonella"])
].copy()

# Select control and salmonella data
data_salm = sccoda_data.mod["coda_salm"]
data_all = sccoda_data.mod["coda"]

# model all three diseases at once
model_all = mod.CompositionalAnalysis(data_all, formula="condition", reference_cell_type="Endocrine")
all_results = model_all.sample_hmc()
all_results.summary()

# Set salmonella infection as "default" category
model_salm_switch_cond = mod.CompositionalAnalysis(
    data_salm,
    formula="C(condition, Treatment('Salmonella'))",
    reference_cell_type="Goblet",
)
switch_results = model_salm_switch_cond.sample_hmc()
switch_results.summary()

# switching reference cell type
model_salm_ref = mod.CompositionalAnalysis(data_salm, formula="condition", reference_cell_type="Enterocyte")
reference_results = model_salm_ref.sample_hmc()
reference_results.set_fdr(est_fdr=0.4)
reference_results.summary()

## result analysis
model_salm = mod.CompositionalAnalysis(data_salm, formula="condition", reference_cell_type="Goblet")
salm_results = model_salm.sample_hmc(num_results=20000)

# extended summary
salm_results.summary_extended(hdi_prob=0.8)

# Run scCODA with each cell type as the reference
cell_types = data_salm.var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["times_credible"]).fillna(0)

for ct in cell_types:
    print(f"Reference: {ct}")

    # Run inference
    model_temp = mod.CompositionalAnalysis(data_salm, formula="condition", reference_cell_type=ct)
    temp_results = model_temp.sample_hmc(num_results=20000)

    # Select credible effects
    cred_eff = temp_results.credible_effects()
    cred_eff.index = cred_eff.index.droplevel(level=0)

    # add up credible effects
    results_cycle["times_credible"] += cred_eff.astype("int")

# Calculate percentages
results_cycle["pct_credible"] = results_cycle["times_credible"] / len(cell_types)
results_cycle["is_credible"] = results_cycle["pct_credible"] > 0.5
print(results_cycle)
