import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 1) Load and prep data
data_path = "./pp.csv"
data = pd.read_csv(data_path)
molecule_counts = data["Molecule_ChEMBL_ID"]\
    .value_counts()\
    .rename("Number of Molecules")
data = data.merge(molecule_counts, left_on="Molecule_ChEMBL_ID", right_index=True)

# Print how many fall in the FDA approved drug-like zone
old_mask = (data["Molecular Weight"] <= 500) & (data["AlogP"].between(1, 4))
print(f"{old_mask.sum():,} / {len(data):,} in FDA-Approved Drug-Like Zone")

# Plot the scatter
plt.figure(figsize=(12, 8))
ax = plt.gca()
sns.scatterplot(
    data=data,
    x="AlogP",
    y="Molecular Weight",
    size="Number of Molecules",
    hue="pChEMBL Value",
    palette="Reds",
    sizes=(10, 300),
    alpha=0.7,
    ax=ax
)

# Define your empirical antimalarial zone
min_logp, max_logp = 2.0, 9.0
min_mw,    max_mw    = 250.0, 530.0
width  = max_logp - min_logp
height = max_mw    - min_mw

# Add a Rectangle patch in data‐coords
antimalarial_rect = mpatches.Rectangle(
    (min_logp, min_mw),   
    width,               
    height,               
    linewidth=0,
    facecolor='green',
    alpha=0.3
)
ax.add_patch(antimalarial_rect)

# Create a custom handle for legend
rect_handle = mpatches.Patch(facecolor='green', alpha=0.3,
                             label='FDA-Approved Drug-Like Zone')

# Reorder legend entries
handles, labels = ax.get_legend_handles_labels()

# Remove any existing entry for label
if 'FDA-Approved Drug-Like Zone' in labels:
    idx = labels.index('FDA-Approved Drug-Like Zone')
    handles.pop(idx)
    labels.pop(idx)

# Append 
handles.append(rect_handle)
labels.append('FDA-Approved Drug-Like Zone')

# Create the legend
ax.legend(
    handles=handles,
    labels=labels,
    title='pIC50 and Assays',
    title_fontproperties={'weight': 'bold'},
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    scatterpoints=1,
    markerscale=0.6,
    labelspacing=1.2
)

# Customize
ax.set_xlabel('AlogP (lipophilicity)', fontsize=12)
ax.set_ylabel('Molecular Weight (Da)', fontsize=12)
ax.set_title('MW vs AlogP, Colored by pIC50 and Sized by # of Molecules', fontsize=14)
ax.grid(True)
plt.tight_layout()

# Save
plt.savefig(
    "./mw_vs_alogp_scat_scaff.png",
    format="png",
    bbox_inches="tight"
)
plt.show()
