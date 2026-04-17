import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.cm import ScalarMappable
import numpy as np

# Global paths
output_dir = "../Analyses/multi-omics/mmvec_outputs/experiment5"
vip_path = "../Data/metabolomics/Run3_10252024/output/VIP_filtered_final_data_method2_11212024_with_shortened_names.csv"
classyfire_path = "../Data/metabolomics/Run3_10252024/info_feature_complete.csv"
# masst_path = "masst_screening_results_05152025_07152025_07192025.csv"
masst_path = "../Data/metabolomics/Run3_10252024/masst_screening_results.csv"

fig_dir_main = "../Figures/Main/Figure_6"
fig_dir_suppl = "../Figures/Supplementary/Suppl_Figure_11"

os.makedirs(fig_dir_main, exist_ok=True)
os.makedirs(fig_dir_suppl, exist_ok=True)

# Comparisons config
comparisons = {
    "Healthy_vs_Acne_L": {
        "title": "MMVEC Microbe–Metabolite Co-Occurrences: H vs Acne L",
        "red_group": "Acne L",
        "blue_group": "Healthy",
        "out_dir": fig_dir_main
    },
    "Healthy_vs_Acne_NL": {
        "title": "MMVEC Microbe–Metabolite Co-Occurrences: H vs Acne NL",
        "red_group": "Acne NL",
        "blue_group": "Healthy",
        "out_dir": fig_dir_suppl
    },
    "Acne_NL_vs_Acne_L": {
        "title": "MMVEC Microbe–Metabolite Co-Occurrences: Acne NL vs Acne L",
        "red_group": "Acne L",
        "blue_group": "Acne NL",
        "out_dir": fig_dir_suppl
    }
}

# Load shared annotation data
VIPs = (
    pd.read_csv(vip_path)
    .drop(columns=["Unnamed: 0.1", "Unnamed: 0"], errors="ignore")
    .set_index("ID")
)
VIPs.index = VIPs.index.astype(str)

compound_name_map = VIPs["Shortened_Compound_Name"].to_dict()

classyfire_info = (
    pd.read_csv(classyfire_path)[["Feature", "ClassyFire#most specific class"]]
    .dropna()
)
classyfire_info["Feature"] = classyfire_info["Feature"].astype(str)
classyfire_map = classyfire_info.set_index("Feature")["ClassyFire#most specific class"].to_dict()

shared_map = {
    compound_name_map[fid]: classyfire_map[fid]
    for fid in compound_name_map if fid in classyfire_map
}

compound_classes = list(shared_map.values())

class_order = [
    'Fatty alcohols', 'Alpha amino acids', 'Alpha amino acids and derivatives',
    'Gamma-glutamyl amino acids', '5-O-methylated flavonoids',
    'N-acyl amines', 'Phenylbenzimidazoles', 'Phenylalanine and derivatives',
    'Organo nitrogen compounds', 'Phenylpropanes',
    'Pyrimidine nucleosides', 'Proline and derivatives',
    'Monosaccharides', '3-O-methylated flavonoids',
    'N-acyl-alpha amino acids', 'Alkyl glycosides',
    'Heteroaromatic compounds', 'Unknown'
]

compound_class_list = sorted(
    set(compound_classes),
    key=lambda x: class_order.index(x) if x in class_order else len(class_order)
)

palette = sns.color_palette("hls", len(compound_class_list))
class_color_map = dict(zip(compound_class_list, palette))
class_colors = [class_color_map[c] for c in compound_classes]

# MASST preprocessing
masst_results = pd.read_csv(masst_path)
masst_results["Category"] = "MASST: " + masst_results["Category"].astype(str)

VIPs_reset = VIPs.reset_index()[["ID", "Shortened_Compound_Name", "RT", "mz"]]
VIPs_reset["Formatted_Name"] = VIPs_reset.apply(
    lambda r: f"{r['Shortened_Compound_Name']} (RT: {r['RT']}, mz: {r['mz']})",
    axis=1
)
compound_id_to_name = dict(zip(VIPs_reset["ID"], VIPs_reset["Formatted_Name"]))

# Plotting function
def plot_mmvec_heatmap(type_key, cfg):

    df_path = f"{output_dir}/{type_key}/conditionals-viz.tsv"
    df = pd.read_csv(df_path, sep="\t").set_index("featureid")

    df.columns = df.columns.str.replace("^g__", "g_", regex=True)

    sorted_cols = [
        'g_Cutibacterium', 'g_Corynebacterium', 'g_Staphylococcus',
        'g_Streptococcus', 'g_Lawsonella', 'g_Micrococcus',
        'g_Lactobacillus', 'g_Kocuria', 'g_Rothia', 'g_Haemophilus'
    ]

    df = df[sorted_cols]
    df = df.loc[df.sum(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(6, 10), dpi=300)

    sns.heatmap(
        df,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        cbar=False,
        linewidths=0.15,
        linecolor="black"
    )

    ax.set_title(cfg["title"], fontsize=22, pad=25)
    ax.tick_params(axis="x", rotation=270, labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_ylabel("")

    out_path = f"{cfg['out_dir']}/heatmap_{type_key}_with_classyfire.png"
    plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out_path}")

# Run all plots
for type_key, cfg in comparisons.items():
    plot_mmvec_heatmap(type_key, cfg)
