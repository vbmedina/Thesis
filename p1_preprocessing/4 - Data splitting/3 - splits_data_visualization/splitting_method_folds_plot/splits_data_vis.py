"""
Data Split Graphic, illustrates:
- 5 splits (K = 1, 2, 3, 4, 5)
- Each split shows 5 equal blocks (each 20% of total data)
- The designated test block shifts from block 5 down to 1 across K fold
- Validation: 10% of the TRAINING set

No requirements
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

# Parametes
TOTAL = 39_624
FOLDS = 5
TRAIN_PCT = 0.70
TEST_PCT = 0.20
VAL_FRAC_OF_TRAIN = 0.10 # 10% of the training set

TEST_BLOCKS = [5, 4, 3, 2, 1]

# Generate a binary matrix where rows=K, cols=blocks 
M = np.zeros((FOLDS, FOLDS), dtype=int)
for k_idx, test_block in enumerate(TEST_BLOCKS):
    M[k_idx, test_block - 1] = 1

# Plot
sns.set_theme(style="white")
fig, ax = plt.subplots(figsize=(14, 3.2))

cmap_reds = sns.color_palette("Reds", as_cmap=True)
test_color = cmap_reds(0.90) # deep red for test
cmap = ListedColormap(["#ffffff", test_color])

sns.heatmap(
    M,
    cmap=cmap,
    cbar=False,
    linewidths=2,
    linecolor="#d1d5db", # grid lines
    ax=ax,
    vmin=0,
    vmax=1,
)

# Off-white background
ax.set_facecolor("#f8fafc")

# Axes labels and ticks
ax.set_xlabel("Blocks (each 20% of dataset)")
ax.set_ylabel("Splits")
ax.set_xticks(np.arange(FOLDS) + 0.5)
ax.set_xticklabels([str(i) for i in range(1,FOLDS + 1)], rotation = 0)
ax.set_yticks(np.arange(FOLDS) + 0.5)
ax.set_yticklabels([f"K = {k}" for k in range(1,FOLDS + 1)], rotation = 0)

# Annotate test cells
for i in range(FOLDS):
    for j in range(FOLDS):
        if M[i, j] == 1:
            ax.text(
                j + 0.5,
                i + 0.5,
                "TEST",
                ha="center",
                va="center",
                fontsize=9,
                color="white",
                fontweight="bold")

# Title and legend numbers
train_n = round(TOTAL * TRAIN_PCT)
test_n = TOTAL - train_n
val_n = round(TOTAL * VAL_FRAC_OF_TRAIN) # 10% of TRAINING set
val_pct_of_total = (val_n / train_n) * 100 # ≈12% of total

title = f"CHEMBL 364 Data Splitting into 5 Folds ({TOTAL:,} total)"
ax.set_title(title, loc="center", fontsize=16, pad=16)

# Legend: indicate validation is carved from training
legend_handles = [
    mpatches.Patch(facecolor=test_color, edgecolor="#cccccc",
                   label=f"Test {int(TEST_PCT*100)}% (≈ {test_n:,})"),
    mpatches.Patch(facecolor="#ffffff", edgecolor="#cccccc",
                   label=f"Train {int(TRAIN_PCT*100)}% (≈ {train_n:,})"),
    mpatches.Patch(facecolor="#ffffff", edgecolor="#666666", hatch="////", linewidth=1.0,
                   label=f"Validation 10% of Total Data (≈ { val_n:,})\n ≈ {val_pct_of_total:.1f}% of Training Data")]

leg = ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1),
                loc="upper left", borderaxespad=0.0)
leg.get_frame().set_facecolor("white")
leg.get_frame().set_edgecolor("#e5e7eb")

plt.tight_layout()
plt.savefig("./p1_preprocessing/4 - Data splitting/3 - splits_data_visualization/splits_data_vis.png", dpi=300, bbox_inches="tight")
plt.show()
