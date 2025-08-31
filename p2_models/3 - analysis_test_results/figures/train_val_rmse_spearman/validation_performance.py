''' Description: Create dot plot with validation results contrasting RMSE to Spearman's rho

Requirements:
    - Non import requirement: Validation data
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
ax.set_title("Validation: Spearman’s rho vs RMSE")

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
            ha="left", va="bottom", fontsize=9)

# Legend outside
ax.legend(title="Split", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)

plt.tight_layout()
plt.savefig("p2_models/3 - analysis_test_results/figures/validation_rmse_spearman/validation_performance.png", dpi=300)
plt.show()
