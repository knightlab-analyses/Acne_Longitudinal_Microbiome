import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib import cm
from matplotlib.lines import Line2D


# ------------------------- Setup ------------------------- #

# Define the color map
colormap = cm.get_cmap("RdBu_r")

# Load conditional probability tables
cond_probs_Healthy_vs_Acne_L = pd.read_csv("../Analyses/multi-omics/mmvec_outputs/experiment5/Healthy_vs_Acne_L/conditionals-viz.tsv", sep='\t', index_col=0)
cond_probs_Healthy_vs_Acne_NL = pd.read_csv("../Analyses/multi-omics/mmvec_outputs/experiment5/Healthy_vs_Acne_NL/conditionals-viz.tsv", sep='\t', index_col=0)
cond_probs_Acne_NL_vs_Acne_L = pd.read_csv("../Analyses/multi-omics/mmvec_outputs/experiment5/Acne_NL_vs_Acne_L/conditionals-viz.tsv", sep='\t', index_col=0)
# cond_probs_AllGroups = pd.read_csv("Outputs/experiment6/conditionals-viz.tsv", sep='\t', index_col=0)

fig_out_main = "../Figures/Main/Figure_6"
fig_out_suppl = "../Figures/Supplementary/Suppl_Figure_11"

# ------------------- Helper Functions -------------------- #

# Function to get top n edges by absolute log_cond_prob (default is 50)
def topN_abs(df, n=50):
    melted = df.reset_index().melt(id_vars='featureid', var_name='microbe', value_name='log_cond_prob')
    return melted.reindex(melted['log_cond_prob'].abs().sort_values(ascending=False).index).head(n)

# Define metabolite node colors based on VIP directionality of association
custom_metabolite_colors_H_vs_Acne_L = {
    '(Iso)leucine': '#dd7966',
    'Tryptophan': '#dd7966',
    'Pyroglutamic\nacid': '#dd7966',
    'Phenylalanine': '#dd7966',
    'C19H36O4\nFatty alcohol': '#dd7966',
    'N-Oleoylethanolamine': '#dd7966',
    'NCGC00380271-01': '#dd7966',
    'Urocanic acid': '#4343a3',
    'Glutamyltyrosine': '#4343a3',
    'Gln-C14:0': '#4343a3'
    # Add more if needed
}

custom_metabolite_colors_H_vs_Acne_NL = {
    '(Iso)leucine': '#6ab2bd',
    'Tryptophan': '#6ab2bd',
    'Pyroglutamic\nacid': '#6ab2bd',
    'Phenylalanine': '#6ab2bd',
    'C19H36O4\nFatty alcohol': '#6ab2bd',
    'N-Oleoylethanolamine': '#6ab2bd',
    'DL-Panthenol': '#4343a3',
    'NCGC00380271-01': '#6ab2bd',
    'Urocanic acid': '#4343a3',
    'Tyrosine': '#6ab2bd'
    # Add more if needed
}

custom_metabolite_colors_Acne_NL_vs_Acne_L = {
    '(Iso)leucine': '#6ab2bd',
    'Tryptophan': '#dd7966',
    'Phenylalanine': '#6ab2bd',
    'N-Oleoylethanolamine': '#dd7966',
    'NCGC00380271-01': '#dd7966',
    'Urocanic acid': '#dd7966',
    'Glutamyltyrosine': '#6ab2bd',
    'Gln-C14:0': '#6ab2bd',
    'Tyrosine': '#6ab2bd',
    'Sinensetin': '#dd7966',
    'xylometazoline': '#dd7966'
    # Add more if needed
}

# Network plotting function
def draw_single_network(top_edges_df, title, outpath, norm_range, label_tweaks=None, layout_seed=42):
    # Clean labels
    top_edges_df['clean_metabolite'] = top_edges_df['featureid'].str.replace(r'\s*\(.*?\)', '', regex=True)
    top_edges_df['clean_metabolite'] = top_edges_df['clean_metabolite'].replace({'leucine': '(Iso)leucine'})
    top_edges_df['microbe_clean'] = top_edges_df['microbe'].str.replace('g__', 'g_', regex=False)

    top_edges_df['clean_metabolite'] = top_edges_df['clean_metabolite'].replace({
        'C19H36O4 Fatty alcohol': 'C19H36O4\nFatty alcohol',
        'Pyroglutamic acid': 'Pyroglutamic\nacid',
        'Sorbitane Monooleate': 'Sorbitane\nMonooleate',
        '3-O-methylated flavonoids': '3-O-methylated\nflavonoids'
    })

    # Create graph
    G = nx.Graph()
    microbes = top_edges_df['microbe_clean'].unique()
    metabolites = top_edges_df['clean_metabolite'].unique()
    G.add_nodes_from(microbes, bipartite=0)
    G.add_nodes_from(metabolites, bipartite=1)

    for _, row in top_edges_df.iterrows():
        G.add_edge(row['microbe_clean'], row['clean_metabolite'], weight=row['log_cond_prob'])

    # Layout using absolute weights
    abs_G = G.copy()
    for u, v, d in abs_G.edges(data=True):
        d['weight'] = abs(d['weight'])  # for layout only

    pos = nx.spring_layout(abs_G, seed=layout_seed, k=0.5)

    # Drawing
    fig, ax = plt.subplots(figsize=(12, 8))


    weights = [d['weight'] for _, _, d in G.edges(data=True)]
    vmax = max(abs(global_min), abs(global_max))
    norm = Normalize(vmin=-vmax, vmax=vmax)

    widths = 2

    nx.draw_networkx_edges(
        G, pos,
        width=widths,
        edge_color=weights,
        edge_cmap=colormap,
        edge_vmin=norm_range[0],
        edge_vmax=norm_range[1],
        alpha=0.8,
        ax=ax
    )

    # Separate nodes
    microbe_nodes = [node for node in G.nodes if node in microbes]
    metabolite_nodes = [node for node in G.nodes if node in metabolites]

    # Node colors
    microbe_colors = ['#ffa505'] * len(microbe_nodes)
    if title == "MMVEC Co-Occurrence Network: Healthy vs Acne Lesional":
        metabolite_colors = [custom_metabolite_colors_H_vs_Acne_L.get(node, '#d3d3d3') for node in metabolite_nodes]
    elif title == "MMVEC Co-Occurrence Network: Healthy vs Acne Non-lesional":
        metabolite_colors = [custom_metabolite_colors_H_vs_Acne_NL.get(node, '#d3d3d3') for node in metabolite_nodes]
    elif title == "MMVEC Co-Occurrence Network: Acne Non-lesional vs Acne Lesional":
        metabolite_colors = [custom_metabolite_colors_Acne_NL_vs_Acne_L.get(node, '#d3d3d3') for node in metabolite_nodes]

    # Draw microbial nodes as squares
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=microbe_nodes,
        node_color=microbe_colors,
        node_size=1000,
        node_shape='s',
        alpha=0.85,
        ax=ax
    )

    # Draw metabolite nodes as circles
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=metabolite_nodes,
        node_color=metabolite_colors,
        node_size=1000,
        node_shape='o',
        alpha=0.85,
        ax=ax
    )


    # Labels
    label_pos = {}
    for node, (x, y) in pos.items():
        if label_tweaks and node in label_tweaks:
            dx, dy = label_tweaks[node]
        else:
            dx, dy = (0.02, 0.02) if node in microbes else (-0.02, -0.02)
        label_pos[node] = (x + dx, y + dy)

    nx.draw_networkx_labels(G, label_pos, font_size=12, font_weight='bold', ax=ax)

    # Title and Colorbar
    ax.set_title(title, fontsize=20)
    ax.axis('off')
    sm = ScalarMappable(norm=norm, cmap=colormap)
    sm.set_array([])

    # Add colorbar in custom position — adjust the values to shrink and center
    cbar_ax = fig.add_axes([0.35, -0.05, 0.3, 0.02])  # [left, bottom, width, height]
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')

    # Label and tick formatting
    cbar.set_label('Log Conditional Probability (Edges)', fontsize=14, fontweight='bold')
    cbar.outline.set_visible(False)

    cbar.ax.tick_params(labelsize=12)

    # ------------------ Node Legend (Bottom Right) ------------------ #
    if title == "MMVEC Co-Occurrence Network: Healthy vs Acne Lesional":
        node_legend_elements = [
            Line2D([0], [0], marker='s', color='w', label='16S taxa',
                markerfacecolor='#ffa505', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Acne L-associated metabolite',
                markerfacecolor='#f16c52', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Healthy-associated metabolite',
                markerfacecolor='#3333B3', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Non-differential metabolite',
                markerfacecolor='#d3d3d3', markersize=16),    
        ]

    elif title == "MMVEC Co-Occurrence Network: Healthy vs Acne Non-lesional":
        node_legend_elements = [
            Line2D([0], [0], marker='s', color='w', label='16S taxa',
                markerfacecolor='#ffa505', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Acne NL-associated metabolite',
                markerfacecolor='#6ab2bd', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Healthy-associated metabolite',
                markerfacecolor='#3333B3', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Non-differential metabolite',
                markerfacecolor='#d3d3d3', markersize=16),    
        ]

    elif title == "MMVEC Co-Occurrence Network: Acne Non-lesional vs Acne Lesional":
        node_legend_elements = [
            Line2D([0], [0], marker='s', color='w', label='16S taxa',
                markerfacecolor='#ffa505', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Acne-NL-associated metabolite',
                markerfacecolor='#6ab2bd', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Acne L-associated metabolite',
                markerfacecolor='#f16c52', markersize=16),
            Line2D([0], [0], marker='o', color='w', label='Non-differential metabolite',
                markerfacecolor='#d3d3d3', markersize=16),    
        ]
    if title == "MMVEC Co-Occurrence Network: Healthy vs Acne Lesional":
        bbox_to_anchor = (0.25, -0.15)
    elif title == "MMVEC Co-Occurrence Network: Healthy vs Acne Non-lesional":
        bbox_to_anchor = (1.0, -0.13)
    else:
        bbox_to_anchor = (1.0, -0.13)

    node_legend = fig.legend(
        handles=node_legend_elements,
        title="Node Key",
        loc='lower right',
        bbox_to_anchor=bbox_to_anchor,
        fontsize=14,
        title_fontsize=14,
        frameon=False
    )
    node_legend.get_title().set_fontweight("bold")

    # Save figure
    plt.tight_layout()
    plt.savefig(outpath, dpi=600, bbox_inches='tight')
    plt.close()

# ------------------- Prepare Top n ------------------- #

# number_of_edges = 48 # top 20%
number_of_edges = 60 # top 25%
# number_of_edges = 72  # top 30%

top_Healthy_vs_Acne_L = topN_abs(cond_probs_Healthy_vs_Acne_L, number_of_edges)
top_Healthy_vs_Acne_NL = topN_abs(cond_probs_Healthy_vs_Acne_NL, number_of_edges)
top_Acne_NL_vs_Acne_L = topN_abs(cond_probs_Acne_NL_vs_Acne_L, number_of_edges)
# top_AllGroups = topN_abs(cond_probs_AllGroups, number_of_edges)

# Compute global normalization range
all_vals = pd.concat([
    top_Healthy_vs_Acne_L['log_cond_prob'],
    top_Healthy_vs_Acne_NL['log_cond_prob'],
    top_Acne_NL_vs_Acne_L['log_cond_prob']
    # top_AllGroups['log_cond_prob']
])
global_min, global_max = all_vals.min(), all_vals.max()



# ------------------- Plot Each Network ------------------- #

label_tweaks = {
    'g_Corynebacterium': (0.0, -0.04),
    'g_Cutibacterium': (0.0, -0.02),
    'g_Staphylococcus': (0.0, -0.02),
    'Sinensetin': (0.0, 0.02),
}

draw_single_network(
    top_Healthy_vs_Acne_L,
    "MMVEC Co-Occurrence Network: Healthy vs Acne Lesional",
    f"{fig_out_main}/Network_Healthy_vs_Acne_L.png",
    (global_min, global_max),
    label_tweaks=label_tweaks
)
 
label_tweaks = {
    'DL-Panthenol': (0.0, 0.03),
    '(Iso)leucine': (0.03, 0.01),
    'Nobiletin': (0.0, 0.01),
    'Sinensetin': (0.0, 0.02),
    'g_Streptococcus': (0.0, -0.01)
}

draw_single_network(
    top_Healthy_vs_Acne_NL,
    "MMVEC Co-Occurrence Network: Healthy vs Acne Non-lesional",
    f"{fig_out_suppl}/Network_Healthy_vs_Acne_NL.png",
    (global_min, global_max),
    label_tweaks=label_tweaks
)

label_tweaks = {
    'g_Corynebacterium': (0.0, -0.02),
    '(Iso)leucine': (0.0, 0.03),
    'g_Cutibacterium': (0.0, -0.02),
    'g_Staphylococcus': (0.0, -0.02)
}

draw_single_network(
    top_Acne_NL_vs_Acne_L,
    "MMVEC Co-Occurrence Network: Acne Non-lesional vs Acne Lesional",
    f"{fig_out_suppl}/Network_Acne_NL_vs_Acne_L.png",
    (global_min, global_max),
    label_tweaks=label_tweaks
)


# draw_single_network(
#     top_AllGroups,
#     "MMVEC Co-Occurrence Network",
#     "Figures/Network_AllGroups.png",
#     (global_min, global_max),
#     label_tweaks=label_tweaks
# )


print("Finished!")

