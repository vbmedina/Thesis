''' Description: Create dot plot with validation results contrasting RMSE to Spearman's rho

Requirements:
    - Non import requirement: Validation data
'''
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

val_rows = [
    {"Split": "Random", "Fold": 1, "RMSE": 0.561, "Spearman": 0.900, "n": 3962},
    {"Split": "Scaffold", "Fold": 1, "RMSE": 0.547, "Spearman": 0.903, "n": 3962},
    {"Split": "Butina", "Fold": 2, "RMSE": 0.553, "Spearman": 0.849, "n": 3962},
    {"Split": "UMAP (k-means)", "Fold": 1, "RMSE": 0.541, "Spearman": 0.900, "n": 3962},
    {"Split": "UMAP (ward)","Fold": 5, "RMSE": 0.566, "Spearman": 0.888, "n": 3962},
]

df = pd.DataFrame(val_rows)
df["Label"] = df["Split"] + " (fold " + df["Fold"].astype(str) + ")"

sns.set_theme(style="whitegrid", context="talk")

palette = sns.color_palette("Reds_r", n_colors=df["Split"].nunique())

plt.figure(figsize=(7.2, 6.2))
ax = sns.scatterplot(
    data=df,
    x="RMSE",
    y="Spearman",
    hue="Split",
    palette=palette,
    s=140,
    edgecolor="black",
    linewidth=0.6)

# Axis labels with padding
LABEL_PAD = 12
ax.set_xlabel("RMSE", labelpad=LABEL_PAD)
ax.set_ylabel("Spearman's rho", labelpad=LABEL_PAD)
ax.set_title("Validation: RMSE against Spearman's rho vs ")
MAJOR_X = 0.01
MAJOR_Y = 0.005
MINOR_Y = 0.0025

xmin, xmax = df["RMSE"].min(), df["RMSE"].max()
ymin, ymax = df["Spearman"].min(), df["Spearman"].max()

def nice_bounds(vmin, vmax, step, pad_steps=1):
    lo = math.floor(vmin/step)*step - pad_steps*step
    hi = math.ceil (vmax/step)*step + pad_steps*step
    return lo, hi

xlo, xhi = nice_bounds(xmin, xmax, MAJOR_X, pad_steps=1)
ylo, yhi = nice_bounds(ymin, ymax, MAJOR_Y, pad_steps=1)

ax.set_xlim(xlo, xhi)
ax.set_ylim(ylo, yhi)

# X: major ticks only, fewer labels
ax.xaxis.set_major_locator(mticker.MultipleLocator(MAJOR_X))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))
ax.xaxis.set_minor_locator(mticker.NullLocator())  # no minor x ticks

# Y: dense ticks (major + minor)
ax.yaxis.set_major_locator(mticker.MultipleLocator(MAJOR_Y))
ax.yaxis.set_minor_locator(mticker.MultipleLocator(MINOR_Y))
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))

# Grid:
ax.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.6)
ax.grid(True, which="minor", linestyle=":",  linewidth=0.5, alpha=0.35)

# Annotate folds
rightmost_idx = df.nlargest(2, "RMSE").index
for idx, r in df.iterrows():
    if idx in rightmost_idx:
        ax.annotate(
            f'fold {int(r["Fold"])}',
            (r["RMSE"], r["Spearman"]),
            xytext=(-8, 6), textcoords="offset points",
            ha="right", va="bottom", fontsize=11)
    else:
        ax.annotate(
            f'fold {int(r["Fold"])}',
            (r["RMSE"], r["Spearman"]),
            xytext=(8, 6), textcoords="offset points",
            ha="left", va="bottom", fontsize=11)

# Legend outside
ax.legend(title="Split", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)

plt.tight_layout()
plt.savefig("p2_models/3 - analysis_test_results/figures/train_val_rmse_spearman/validation_performance.png", dpi=300)
plt.show()