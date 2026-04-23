"""
Sebum levels: Acne vs. Healthy.

One value per subject (mean across their samples), grouped by cohort.
Comparison via Mann-Whitney U (non-parametric, robust to small n and
outliers — consistent with the Spearman choice in the age-vs-sebum plot).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

INFILE = "../Metadata/metadata_final_22102024_with_age.tsv"
OUTFILE = "../Figures/Supplementary/Suppl_Figure_13/sebum_by_cohort.png"

df = pd.read_csv(INFILE, sep="\t")
dat = df[["subject_ID", "sebumeter_casual_level_mean_forehead_left", "cohort"]].dropna()
dat = dat.rename(columns={"sebumeter_casual_level_mean_forehead_left": "sebum"})

per_subj = (
    dat.groupby("subject_ID")
    .agg(sebum=("sebum", "mean"), cohort=("cohort", "first"))
    .reset_index()
)

acne = per_subj.loc[per_subj["cohort"] == "acne", "sebum"].values
healthy = per_subj.loc[per_subj["cohort"] == "control", "sebum"].values

# Mann-Whitney U (two-sided) — appropriate for small, non-normal samples
u_stat, p_mw = stats.mannwhitneyu(acne, healthy, alternative="two-sided")

print(f"Healthy (n = {len(healthy)}): median = {np.median(healthy):.1f}, "
      f"mean = {healthy.mean():.1f}")
print(f"Acne    (n = {len(acne)}): median = {np.median(acne):.1f}, "
      f"mean = {acne.mean():.1f}")
print(f"Mann-Whitney U = {u_stat:.1f}, p = {p_mw:.4g}")

# --- Figure ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 6))

groups = ["control", "acne"]
labels = ["Healthy", "Acne"]
colors = {"control": "#2c7fb8", "acne": "#FA8072"}
positions = [1, 2]

# Boxplots (no outlier markers — we overlay all points)
bp = ax.boxplot(
    [healthy, acne],
    positions=positions,
    widths=0.55,
    patch_artist=True,
    showfliers=False,
    medianprops=dict(color="black", linewidth=2),
    boxprops=dict(linewidth=1.2),
    whiskerprops=dict(linewidth=1.2),
    capprops=dict(linewidth=1.2),
)
for patch, g in zip(bp["boxes"], groups):
    patch.set_facecolor(colors[g])
    patch.set_alpha(0.35)
    patch.set_edgecolor(colors[g])

# Jittered points
rng = np.random.default_rng(0)
for pos, g, values in zip(positions, groups, [healthy, acne]):
    jitter = rng.uniform(-0.09, 0.09, size=len(values))
    ax.scatter(np.full(len(values), pos) + jitter, values,
               color=colors[g], s=70, zorder=3,
               edgecolor="white", linewidth=0.8)

# Significance bracket
y_max = per_subj["sebum"].max()
y_min = per_subj["sebum"].min()
y_span = y_max - y_min
bar_y = y_max + y_span * 0.08
tick = y_span * 0.02
ax.plot([1, 1, 2, 2], [bar_y, bar_y + tick, bar_y + tick, bar_y],
        color="black", linewidth=1.2)
p_label = f"p = {p_mw:.3f}" if p_mw >= 0.001 else "p < 0.001"
ax.text(1.5, bar_y + tick * 1.5, f"Mann–Whitney  {p_label}",
        ha="center", va="bottom", fontsize=11)

ax.set_xticks(positions)
ax.set_xticklabels([f"Healthy\n(n = {len(healthy)})", f"Acne\n(n = {len(acne)})"], fontsize=13)
ax.set_ylabel("Sebum Level (μg/cm²)", fontsize=14)
ax.tick_params(axis="y", labelsize=12)
ax.set_title("Sebum Levels by Subject Mean: Acne vs. Healthy", fontsize=14, pad=12)
ax.set_ylim(y_min - y_span * 0.08, bar_y + tick * 6)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(True, axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(OUTFILE, dpi=300, bbox_inches="tight")
print(f"\nSaved: {OUTFILE}")
