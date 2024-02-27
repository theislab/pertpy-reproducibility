import scanpy as sc

# plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica"]

pt_blue = "#34669A"
pt_red = "#D2455E"
pt_orange = "#F37800"
pt_gold = "#FFD700"


def figure_journal_basic():
    sc.set_figure_params(
        dpi=300,
        dpi_save=300,
        figsize=(4, 4),
        frameon=False,
        facecolor=None,
        transparent=True,
        vector_friendly=True,
    )


def figure_jupyter():
    sc.set_figure_params(
        dpi=80,
        dpi_save=150,
        figsize=(8, 8),
        frameon=False,
        facecolor="white",
        transparent=False,
        vector_friendly=True,
    )


mhv68_pcls_ct_color_map = {}

mcfarland_ct_color_map = {}

norman_ct_color_map = {}

norman_gene_program_color_map = {
    "Control": pt_gold,
    "Erythroid": "Green",
    "G1 cell cycle": pt_blue,
    "G1 cell cycle (KNN imputed)": "#648fbd",
    "Granulocyte apoptosis": "LightSeaGreen",
    "Megakaryocyte": "Purple",
    "Pioneer factors": pt_red,
    "Pioneer factors (KNN imputed)": "#e39aa8",
    "Pro-growth": pt_orange,
    "Pro-growth (KNN imputed)": "#e8ac72",
    "Unknown": "LightGrey",
}
