"""
Correlation between age and sebumeter (casual level, forehead left).

The input has 266 rows but only 20 unique subjects (each subject
contributes multiple samples). We report correlations both per-row
and per-subject. The per-subject estimate is the statistically cleaner
one because it does not treat pseudo-replicates as independent.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

INFILE = "../Metadata/metadata_final_22102024_with_age.tsv"
OUTFILE = "../Figures/Supplementary/Suppl_Figure_13/age_sebum_correlation.png"

df = pd.read_csv(INFILE, sep="\t")
dat = df[["subject_ID", "age", "sebumeter_casual_level_mean_forehead_left", "cohort"]].dropna()
dat = dat.rename(columns={"sebumeter_casual_level_mean_forehead_left": "sebum"})

# --- Per-row correlation (n = 266 samples) --------------------------------
r_row, p_row = stats.pearsonr(dat["age"], dat["sebum"])
rho_row, p_rho_row = stats.spearmanr(dat["age"], dat["sebum"])

# --- Per-subject correlation (n = 20 subjects) ----------------------------
per_subj = (
    dat.groupby("subject_ID")
    .agg(age=("age", "first"),
         sebum=("sebum", "mean"),
         cohort=("cohort", "first"))
    .reset_index()
)
r_sub, p_sub = stats.pearsonr(per_subj["age"], per_subj["sebum"])
rho_sub, p_rho_sub = stats.spearmanr(per_subj["age"], per_subj["sebum"])

print(f"Per-row (n = {len(dat)}):")
print(f"  Pearson  r   = {r_row:.3f}, p = {p_row:.4g}")
print(f"  Spearman rho = {rho_row:.3f}, p = {p_rho_row:.4g}\n")
print(f"Per-subject (n = {len(per_subj)}):")
print(f"  Pearson  r   = {r_sub:.3f}, p = {p_sub:.4g}")
print(f"  Spearman rho = {rho_sub:.3f}, p = {p_rho_sub:.4g}")

# --- Figure ---------------------------------------------------------------
# Grey points: individual samples. Blue points: per-subject means.
# Regression line + 95% CI is fit on the per-subject means (correct
# level of analysis here).

fig, ax = plt.subplots(figsize=(7, 7))

# Regression line + 95% CI on per-subject means
x = per_subj["age"].values.astype(float)
y = per_subj["sebum"].values.astype(float)
slope, intercept = np.polyfit(x, y, 1)
x_line = np.linspace(x.min() - 0.5, x.max() + 0.5, 100)
y_line = slope * x_line + intercept

# 95% confidence band around the regression line
n = len(x)
y_hat = slope * x + intercept
resid = y - y_hat
s_err = np.sqrt(np.sum(resid ** 2) / (n - 2))
x_mean = x.mean()
t_val = stats.t.ppf(0.975, n - 2)
ci = t_val * s_err * np.sqrt(
    1.0 / n + (x_line - x_mean) ** 2 / np.sum((x - x_mean) ** 2)
)
ax.fill_between(x_line, y_line - ci, y_line + ci, color="grey", alpha=0.25)
ax.plot(x_line, y_line, color="black", lw=2)

# Subject means, coloured by cohort
colors = {"acne": "#FA8072", "control": "#2c7fb8"}  # salmon red / blue
labels = {"acne": "Acne", "control": "Healthy"}
for coh, grp in per_subj.groupby("cohort"):
    ax.scatter(grp["age"], grp["sebum"], color=colors[coh], s=70,
               zorder=3, label=labels[coh], edgecolor="white", linewidth=0.8)

annot = (
    f"Spearman ρ = {rho_sub:.2f}\n"
    f"p = {p_rho_sub:.3f}\n"
    f"n = {len(per_subj)}"
)
ax.text(0.98, 0.97, annot, transform=ax.transAxes,
        ha="right", va="top", fontsize=12,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="lightgrey", alpha=0.9))

ax.set_xticks(range(int(dat["age"].min()), int(dat["age"].max()) + 1))
ax.set_xlabel("Age (years)", fontsize=16)
ax.set_ylabel("Sebum Level (μg/cm²)", fontsize=16)
ax.tick_params(axis="both", labelsize=14)
ax.set_title("Correlation of Sebum Levels by Subject Mean with Age", fontsize=16, loc="center", pad=12)
ax.legend(loc="lower right", frameon=False, fontsize=13)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTFILE, dpi=300, bbox_inches="tight")
print(f"\nSaved: {OUTFILE}")
