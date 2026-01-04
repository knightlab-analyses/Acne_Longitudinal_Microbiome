import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.cm import ScalarMappable
import matplotlib.patches as patches
import numpy as np

# --------------------------- Setup --------------------------- #
output_dir = "Outputs/experiment5"
vip_path = "VIP_filtered_with_shortened_names.tsv"
classyfire_path = "../Data/metabolomics/Run3_10252024/info_feature_complete.csv"
masst_path = "masst_screening_results_05152025_07152025_07192025.csv"
fig_dir = "Figures"

## Uncomment to run the following:
# type_key = "Healthy_vs_Acne_L"
# title = "MMVEC Microbe-Metabolite Co-Occurrences: H vs Acne L"
# red_group = "Acne L"
# blue_group = "Healthy"

# type_key = "Healthy_vs_Acne_NL"
# title = "MMVEC Microbe-Metabolite Co-Occurrences: H vs Acne NL"
# red_group = "Acne NL"
# blue_group = "Healthy"

type_key = "Acne_NL_vs_Acne_L"
title = "MMVEC Microbe-Metabolite Co-Occurrences: Acne NL vs Acne L"
red_group = "Acne L"
blue_group = "Acne NL"

os.makedirs(fig_dir, exist_ok=True)

# ----------------------- Load Data --------------------------- #
VIPs = pd.read_csv(vip_path, sep='\t').drop(columns=["Unnamed: 0.1", "Unnamed: 0"], errors='ignore')
VIPs.set_index("ID", inplace=True)
VIPs.index = VIPs.index.astype(str)
compound_name_map = VIPs["Shortened_Compound_Name"].to_dict()

# Load ClassyFire annotations
classyfire_info = pd.read_csv(classyfire_path)[["Feature", "ClassyFire#most specific class"]].dropna()
classyfire_info["Feature"] = classyfire_info["Feature"].astype(str)
classyfire_map = classyfire_info.set_index("Feature")["ClassyFire#most specific class"].to_dict()

# Merge VIP and ClassyFire
shared_map = {compound_name_map[fid]: classyfire_map[fid] for fid in compound_name_map if fid in classyfire_map}
compound_class = list(shared_map.values())

class_order = [
    'Fatty alcohols', 'Alpha amino acids', 'Alpha amino acids and derivatives', 'Gamma-glutamyl amino acids',
    '5-O-methylated flavonoids', 'N-acyl amines', 'Phenylbenzimidazoles', 'Phenylalanine and derivatives',
    'Organo nitrogen compounds', 'Phenylpropanes', 'Pyrimidine nucleosides', 'Proline and derivatives',
    'Monosaccharides', '3-O-methylated flavonoids', 'N-acyl-alpha amino acids', 'Alkyl glycosides',
    'Heteroaromatic compounds', 'Unknown'
]
compound_class_list = sorted(set(compound_class), key=lambda x: class_order.index(x) if x in class_order else len(class_order))
palette = sns.color_palette("hls", len(compound_class_list))
class_color_map = dict(zip(compound_class_list, palette))
class_colors = [class_color_map[c] for c in compound_class]

# ----------------------- Main Heatmap ------------------------ #
df_path = f"{output_dir}/{type_key}/conditionals-viz.tsv"
df = pd.read_csv(df_path, sep='\t').set_index("featureid")
sorted_rows = df.sum(axis=1).sort_values(ascending=False).index.tolist()
sorted_cols = ['g_Cutibacterium', 'g_Corynebacterium', 'g_Staphylococcus',
               'g_Streptococcus', 'g_Lawsonella', 'g_Micrococcus',
               'g_Lactobacillus', 'g_Kocuria', 'g_Rothia', 'g_Haemophilus']
df.columns = df.columns.str.replace("^g__", "g_", regex=True)
df = df.reindex(sorted_rows).dropna(how="all")
df = df[sorted_cols]


# ------------------ Setup Plotting -------------------------- #
fig, ax = plt.subplots(figsize=(6, 10), dpi=300)
vmin, vmax = -3, 3
cmap = "RdBu_r"

sns.heatmap(df, ax=ax, cmap=cmap, center=0, vmin=vmin, vmax=vmax,
            annot=False, xticklabels=True, yticklabels=True, cbar=False, linewidths=0.15, linecolor='black')

ax.set_ylabel("")
ax.set_title(title, fontsize=22, pad=25)
ax.tick_params(axis="x", rotation=270, labelsize=14)
ax.tick_params(axis="y", rotation=0, labelsize=14)
for spine in ax.spines.values():
    spine.set_visible(True)

# ------------------- Metabolite Direction ------------------- #
direction_colors = {
    'Tyrosine (RT: 0.6720225, mz: 182.0809391)': '#000000', 'Pyroglutamic acid (RT: 0.64796996, mz: 130.0497016)': '#f16c52', 
    'DL-Panthenol (RT: 1.1194918, mz: 206.1384472)': '#000000', 'Phenylalanine (RT: 0.9777797, mz: 166.086109)': '#f16c52', 
    '(Iso)leucine (RT: 0.71511286, mz: 132.1017009)': '#f16c52', 'C19H36O4 Fatty alcohol (RT: 7.4214296, mz: 311.2573967)': '#f16c52',
    'Tryptophan (RT: 2.0760007, mz: 205.0969182)': '#f16c52', 'Cyclo(his-pro)': '#000000', 'Phenylbenzimidazole sulfonic acid': '#f16c52',
    'NCGC00380271-01 (RT: 3.5891602, mz: 505.2636006)': '#f16c52', 'Uridine (RT: 0.65536463, mz: 245.0765944)': '#000000', 
    'Sorbitane Monooleate (RT: 6.1563983, mz: 428.3728208)': '#000000', 'N-Oleoylethanolamine (RT: 7.1323557, mz: 326.3050565)': '#f16c52', 
    'Nobiletin (RT: 5.258404, mz: 403.138254)': '#000000', '3-O-methylated flavonoids (RT: 5.410203, mz: 433.1489336)': '#000000',
    'xylometazoline (RT: 4.3418465, mz: 245.2008949)': '#000000', 'Urocanic acid (RT: 0.6406374, mz: 139.0410336)': '#3333B3',
    'Sinensetin (RT: 5.0247345, mz: 373.127569)': '#000000', 'Glutamyltyrosine (RT: 1.2036344, mz: 311.1238575)': '#3333B3', 
    'Gln-C14:0 (RT: 5.845114, mz: 357.2754847)': '#3333B3'
}

for label in ax.get_yticklabels():
    metabolite = label.get_text()
    label.set_color(direction_colors.get(metabolite, "black"))

# --------------------- MASST Strip (with proper alignment) -------------------------- #
masst_results = pd.read_csv(masst_path)
masst_results['Category'] = 'MASST: ' + masst_results['Category'].astype(str)
unique_compounds = VIPs.index.unique()
unique_categories = masst_results['Category'].dropna().unique()
binary_df = pd.DataFrame(0, index=unique_compounds, columns=unique_categories)
for _, row in masst_results.iterrows():
    binary_df.loc[str(row['Compound_ID']), row['Category']] = 1

# Create formatted name using VIPs info to match df.index
VIPs_reset = VIPs.reset_index()[['ID', 'Shortened_Compound_Name', 'RT', 'mz']]
VIPs_reset['Formatted_Name'] = VIPs_reset.apply(
    lambda row: f"{row['Shortened_Compound_Name']} (RT: {row['RT']}, mz: {row['mz']})", axis=1
)
compound_id_to_formatted_name = dict(zip(VIPs_reset['ID'], VIPs_reset['Formatted_Name']))
binary_df['Formatted_Name'] = binary_df.index.map(compound_id_to_formatted_name)
binary_df = binary_df.dropna(subset=['Formatted_Name'])
binary_df.set_index('Formatted_Name', inplace=True)

# Align with heatmap
binary_df = binary_df.loc[binary_df.index.intersection(df.index)]
binary_df = binary_df.reindex(df.index)

binary_df.to_csv("binary_df.csv")

category_colors = {
    'MASST: microbes': '#76325F66',
    'MASST: microbiome': '#F0897666',
    'MASST: personal care product': '#E8BAD866',
    'MASST: tissue': '#393A9266',
    'MASST: plants': '#45855066',
    'MASST: food': '#3391C166'
}
strip_width = 0.08

for i, category in enumerate(binary_df.columns):
    binary = binary_df[category].fillna(0).astype(int).values
    left = 1.1 + (i * strip_width)
    cmap_bin = ListedColormap(["#ffffff", category_colors.get(category, "#cccccc66")])

    masst_ax = ax.inset_axes([left, 0, strip_width, 1], transform=ax.transAxes)
    
    masst_ax.imshow(
        np.array(binary).reshape(-1, 1),
        aspect='auto',
        extent=[0, 1, 0, len(df)],
        cmap=cmap_bin,
        vmin=0,
        vmax=1,
        origin='lower'
                )
    
    
    for y in range(1, len(df)):
        masst_ax.hlines(y, xmin=0, xmax=1, colors='black', linewidth=0.5)

    masst_ax.set_ylim(len(df), 0)
    masst_ax.set_xticks([0.5])
    masst_ax.set_xticklabels([category], rotation=270, ha='center', fontsize=14)
    masst_ax.set_yticks([])
    masst_ax.tick_params(length=0)

    for spine in masst_ax.spines.values():
        spine.set_visible(True)

# ------------------ ClassyFire Strip ------------------------ #
color_indices = np.arange(len(df)).reshape(-1, 1)
class_cmap = ListedColormap(class_colors)

classy_ax = ax.inset_axes([1.01, 0, 0.08, 1], transform=ax.transAxes)
classy_ax.imshow(
    color_indices,
    aspect='auto',
    cmap=class_cmap,
    extent=[0, 1, 0, len(df)],
    vmin=0,
    vmax=len(class_colors) - 1
)
for y in range(1, len(df)):
    classy_ax.hlines(y, xmin=0, xmax=1, colors='black', linewidth=0.5)
classy_ax.set_ylim(len(df), 0)
classy_ax.set_xticks([])
classy_ax.set_yticks([])
classy_ax.set_xlabel("ClassyFire MSC", fontsize=14, rotation=270, ha="center")
classy_ax.xaxis.set_label_position("bottom")
for spine in classy_ax.spines.values():
    spine.set_visible(True)

# ------------------------ Legends --------------------------- #
norm = Normalize(vmin=vmin, vmax=vmax)
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar_ax = fig.add_axes([-0.8, -0.29, 0.75, 0.03])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
cbar.set_label("Log Conditional Probability", fontsize=14, labelpad=10, fontweight="bold", ha="center")
cbar.ax.tick_params(axis="x", labelsize=10, direction="out", rotation=0)
cbar.outline.set_visible(False)

legend_patches = [
    plt.Line2D([0], [0], marker="s", color="w", label=label,
               markerfacecolor=color, markersize=10)
    for label, color in class_color_map.items()
]
legend = fig.legend(
    handles=legend_patches,
    title="ClassyFire Most Specific Class",
    bbox_to_anchor=(1.25, -0.15),
    fontsize=12,
    ncol=2,
    title_fontsize=14,
    frameon=False
)
legend.get_title().set_fontweight("bold")

fig.text(-0.1, 0.025, "Metabolite Direction", fontsize=14, ha="right", fontweight="bold")
fig.text(-0.1, 0, f"Red = {red_group}-associated", color="#f16c52", fontsize=14, ha="right")
fig.text(-0.1, -0.025, f"Blue = {blue_group}-associated", color="#3333B3", fontsize=14, ha="right")
fig.text(-0.1, -0.05, f"Black = Not differential", color="black", fontsize=14, ha="right")

# ---------------------- Save Plot ---------------------------- #
fig_path = f"{fig_dir}/heatmap_{type_key}_with_classyfire.png"
plt.savefig(fig_path, dpi=600, bbox_inches="tight")
plt.close()
print(f"✅ Saved figure: {fig_path}")
