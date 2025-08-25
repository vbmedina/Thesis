"""
Chemprop (D-MPNN) Early Stopping Visualization

In Chemprop's D-MPNN, early stopping watches the validation metric. It keeps the best-so-far checkpoint (lowest validation 
loss or highest val score). If the validation metric doesn't improve in patience epochs, training stops, and the model 
restores the best checkpoint. Then it will evaluate the test set once with the restored model. (1)

What this figure shows:
- Loss curve on synthetic data.
- Training loss (gray, dashed) usually keeps decreasing.
- Validation loss (red) decreases then rises (overfitting).
- Best epoch (minimum validation loss) is marked with a red dot and dotted red vline.
- Training stops later when patience expires.
- The red shaded region between is the no-improvement window (patience).
- The test metric is evaluated once at the best epoch (yellow star).

References (code/data in paper):
1) CHEMPROP: https://github.com/chemprop/chemprop
"""

import numpy as np
import matplotlib.pyplot as plt

# Synthetic demo curves 
rng = np.random.default_rng(7)
epochs = np.arange(1, 151)

# Training loss: decay + noise
train_loss = 1.5 * np.exp(-epochs / 40) + 0.02 * rng.normal(size=epochs.size) + 0.10

# Validation loss: overfitting curve + noise
val_loss = 1.6 * np.exp(-epochs / 45) + 0.02 * rng.normal(size=epochs.size) + 0.12
val_loss += 0.0006 * np.maximum(0, (epochs - 85)) ** 1.6

best_idx = int(np.argmin(val_loss))
best_epoch = int(epochs[best_idx])

# Patience setting
patience = 15
stop_epoch = min(epochs[-1], best_epoch + patience)

# A test measurement at the best epoch
test_at_best = val_loss[best_idx] + 0.03

# Plotting parameters colors
train_color = "0.35"
val_color   = "#b91c1c"
star_color  = "#facc15"
shade_color = "#fecaca"
stop_color  = "#b91c1c"

plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "-"})

fig, ax = plt.subplots(figsize=(10, 5.2))

# Curves
ax.plot(epochs, train_loss, linestyle="--", linewidth=2.0,
        color=train_color, label="Train loss")
ax.plot(epochs, val_loss, linewidth=2.6, color=val_color,
        marker="o", markevery=6, markersize=4, label="Validation loss")

# Best epoch markers
ax.scatter([best_epoch], [val_loss[best_idx]], s=90, facecolor=val_color,
           edgecolors="white", linewidth=0.8, zorder=5, label="Best (min validation)")
ax.axvline(best_epoch, color=val_color, linestyle=":", linewidth=1.8)

# No-improvement window (patience) and stop epoch
ax.axvspan(best_epoch, stop_epoch, color=shade_color, alpha=0.25,
           label=f"No improvement for {patience} epochs")
ax.axvline(stop_epoch, color=stop_color, linestyle="--", linewidth=1.8)

# Training stops annotation
ax.text(
    stop_epoch + 0.6,
    0.50,
    "training stops",
    transform=ax.get_xaxis_transform(),  # x in data coords, y in axes fraction
    rotation=90,
    va="center", ha="left",
    color=stop_color, fontsize=9,
    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.9, edgecolor="none"),
    clip_on=False, zorder=6,
)

# Test at best epoch (single yellow star)
ax.scatter([best_epoch], [test_at_best], s=200, marker="*",
           facecolor=star_color, edgecolors="white", linewidth=0.9,
           zorder=6, label="Test (evaluated at best epoch)")

# Labels
ax.set_xlabel("Epoch")
ax.set_ylabel("Loss")
ax.set_title("Chemprop's D-MPNN: Early Stopping Driven by Validation (Synthetic Data)", loc="center", fontsize=16, pad=16)

# Legend
leg = ax.legend(frameon=True, framealpha=0.96, facecolor="white", loc="upper right")
leg.get_frame().set_edgecolor("#e5e7eb")

# Final layout and save
fig.tight_layout()
plt.savefig("./p1_preprocessing/5 - Data visualizations after splits /loss_dmpnn_early_stopping_visual/chemprop_dmpnn_early_stopping_validation_driven.png", dpi=300, bbox_inches="tight")
plt.show()
