''' Description: This script is used to generate a histogram of pIC50 distributions with activity thresholds, highlighting 
inactive as >6.0, active as >=6.0 and <7.5, and high-potency as >=7.5. These activity thresholds align with established 
ChEMBL-based antimalarial classification schemes, where pIC50 >= 6.0 represents the minimum threshold for biological activity 
(1). 

References: 
1. pIC50 6.0 minimum threshold for biological activity: https://ieeexplore.ieee.org/abstract/document/10469118?casa_token=McAUgdGGNr0AAAAA:lJBsw0r76ODO4vqL-euLTiJYpqNVC8SbmQSOSbhjZZfGM8ExWeY84Ipoieem6uU2lwb9vN6b

Requirements:
1) final_data.csv - from Step 2 if pipeline.
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors, patches, lines

# Path
csv_path = "./p1_preprocessing/3 - Data visualizations before splits/pIC50_distribution_histogram/final_data_copy.csv"

# Thresholds
main_threshold = 6.0     
upper_threshold = 7.5    
bins = 60
figsize = (7, 4)

# DF
df = pd.read_csv(csv_path, low_memory=False)

# Category Counts
inactive_count = (df["pIC50"] < main_threshold).sum()
active_count = ((df["pIC50"] >= main_threshold) & (df["pIC50"] < upper_threshold)).sum()
high_count = (df["pIC50"] >= upper_threshold).sum()

# Binning
counts, bin_edges = np.histogram(df["pIC50"].dropna(), bins=bins)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Colour Mapping
norm = colors.Normalize(vmin=counts.min(), vmax=counts.max())
cmap = cm.get_cmap("Reds")
bar_colors = cmap(norm(counts))

# Plot
fig, ax = plt.subplots(figsize=figsize)
ax.bar(bin_centers, counts, width=np.diff(bin_edges), align="center",
       edgecolor="black", color=bar_colors)

# Threshold Lines
lt_active = ax.axvline(main_threshold,  ls="-", lw=2.0, color="black")
lt_high   = ax.axvline(upper_threshold, ls="-.", lw=2.0, color="black")

# Legend Patches
patch_inactive = patches.Patch(color=cmap(norm(counts.max() * 0.9)),
                               label=f"Inactive (<{main_threshold}): {inactive_count:,}")
patch_active   = patches.Patch(color=cmap(norm(counts.max() * 0.6)),
                               label=f"Active (≥{main_threshold} & <{upper_threshold}): {active_count:,}")
patch_high     = patches.Patch(color=cmap(norm(counts.max() * 0.3)),
                               label=f"High-potency (≥{upper_threshold}): {high_count:,}")

# legend Lines
line_active = lines.Line2D([], [], color="black", lw=2, linestyle="-",
                           label=f"Active threshold (pIC50={main_threshold})")
line_high   = lines.Line2D([], [], color="black", lw=2, linestyle="-.",
                           label=f"High-potency threshold (pIC50={upper_threshold})")

# Handles
handles = [patch_inactive, patch_active, patch_high, line_active, line_high]
ax.legend(handles=handles, frameon=True, fontsize=9, title="Classes & thresholds")

# Customizations
ax.set_xlabel("pIC50")
ax.set_ylabel("Count")
ax.set_title("pIC50 distribution with activity thresholds")

# Save and Show
plt.tight_layout()
plt.savefig("./p1_preprocessing/3 - Data visualizations before splits/pIC50_distribution_histogram/pIC50_histogram.png", dpi=300)

# Print summary
print(f"Inactive: {inactive_count:,}, Active: {active_count:,}, High-potency: {high_count:,}")
print(f"Active threshold: {main_threshold}, High-potency threshold: {upper_threshold}")
