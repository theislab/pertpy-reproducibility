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
