# import matplotlib.pyplot as plt
# import numpy as np

# # Data for each preprocessing stage in desired order (raw at top)
# stages = [
#     'Raw CSV',
#     'Stage 1',
#     'Stage 2',
#     'Stage 3',
#     'Stage 4',
#     'Stage 5',
#     'Final Dataset'
# ]
# counts = [
#     96923,
#     96923 * 0.79,
#     96923 * 0.65,
#     96923 * 0.57,
#     96923 * 0.57,
#     43482,
#     43482
# ]

# # Generate a blue gradient palette
# cmap = plt.get_cmap('Blues')
# colors = cmap(np.linspace(0.4, 0.8, len(stages)))

# # Create a horizontal bar chart
# fig, ax = plt.subplots(figsize=(6, 8))
# bars = ax.barh(stages, counts, color=colors, edgecolor='black')

# # Annotate counts on bars
# for bar, count in zip(bars, counts):
#     ax.text(count + 2000, bar.get_y() + bar.get_height() / 2,
#             f'{int(count):,}', va='center', fontsize=9)

# # Styling
# ax.set_xlabel('Number of entries', fontsize=12)
# ax.set_title('Dataset Size After Each Preprocessing Stage', fontsize=14, pad=15)

# # Remove spines for a clean look
# for spine in ['top', 'right']:
#     ax.spines[spine].set_visible(False)

# # Remove ticks
# ax.xaxis.set_ticks_position('none')
# ax.yaxis.set_ticks_position('none')

# # Ensure Raw CSV is at the top
# ax.invert_yaxis()

# plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/prepros_viz.svg",
#     format="svg",
#     bbox_inches="tight")
# plt.show()    

import matplotlib.pyplot as plt
import numpy as np

# Data for each preprocessing stage in desired order (raw at top)
stages = [
    'Raw CSV',
    'Stage 1',
    'Stage 2',
    'Stage 3',
    'Stage 4',
    'Stage 5',
    'Final Dataset'
]
counts = [
    96923,
    96923 * 0.79,
    96923 * 0.65,
    96923 * 0.57,
    96923 * 0.57,
    43482,
    43482
]

# Generate a blue gradient palette
cmap = plt.get_cmap('Blues')
colors = cmap(np.linspace(0.4, 0.8, len(stages)))

# Create a wider, shorter horizontal bar chart
fig, ax = plt.subplots(figsize=(10, 4))
bars = ax.barh(stages, counts, color=colors, edgecolor='none')

# Annotate counts on bars
for bar, count in zip(bars, counts):
    ax.text(count + 2000, bar.get_y() + bar.get_height() / 2,
            f'{int(count):,}', va='center', fontsize=9)

# Styling
ax.set_xlabel('Number of entries', fontsize=12)
ax.set_title('Dataset Size After Each Preprocessing Stage', fontsize=14, pad=10)
ax.grid(axis='x', linestyle='--', alpha=0.5)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Invert y-axis so Raw CSV is at top
ax.invert_yaxis()

plt.savefig("/Users/victoriamedina/Thesis_Project/Thesis/Visualizations/prepros_viz.png",
    format="png",
    bbox_inches="tight", dpi=300)
plt.show()
