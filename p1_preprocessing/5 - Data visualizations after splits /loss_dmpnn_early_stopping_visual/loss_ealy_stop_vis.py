"""
Chemprop (D-MPNN) Early Stopping Visualization

In Chemprop/D-MPNN, early stopping watches the validation metric. It keeps the best-so-far checkpoint (lowest validation loss 
or highest val score). If the validation metric doesn't improve in patience epochs, training stops, and the model restores 
the best checkpoint. Then it will evaluate the test set once with the restored model.

What this figure shows:
- Training loss (gray, dashed) usually keeps decreasing.
- Validation loss (red) decreases then rises (overfitting).
- Best epoch (minimum validation loss) is marked with a red dot and dotted red vline.
- Training stops later when patience expires.
- The red shaded region between is the no-improvement window (patience).
- The test metric is evaluated once at the best epoch (yellow star).

Save: PNG (uncomment SVG for vector).
"""


import numpy as np
import matplotlib.pyplot as plt

# --------------------------
# 1) Synthetic demo curves (replace with real logs if available)
# --------------------------
rng = np.random.default_rng(7)
epochs = np.arange(1, 151)

# Training loss: monotonic-ish decay + small noise
train_loss = 1.5 * np.exp(-epochs / 40) + 0.02 * rng.normal(size=epochs.size) + 0.10

# Validation loss: decays then rises (overfitting)
val_loss = 1.6 * np.exp(-epochs / 45) + 0.02 * rng.normal(size=epochs.size) + 0.12
val_loss += 0.0006 * np.maximum(0, (epochs - 85)) ** 1.6  # upturn after ~85

best_idx = int(np.argmin(val_loss))
best_epoch = int(epochs[best_idx])

# Patience setting (Chemprop-style)
patience = 15
stop_epoch = min(epochs[-1], best_epoch + patience)

# A single test measurement at the best epoch (demo value near val best)
test_at_best = val_loss[best_idx] + 0.03

# --------------------------
# 2) Plot
# --------------------------
# Colors
train_color = "0.35"          # gray
val_color   = "#b91c1c"       # red
star_color  = "#facc15"       # yellow
shade_color = "#fecaca"       # light red for no-improvement window
stop_color  = "#6b7280"       # gray for stop vline

plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "-"
})

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
ax.text(stop_epoch + 1, ax.get_ylim()[0] + 0.03*(ax.get_ylim()[1]-ax.get_ylim()[0]),
        "training stops", color=stop_color, fontsize=9, rotation=90, va="center")


# Test at best epoch (single yellow star)
ax.scatter([best_epoch], [test_at_best], s=200, marker="*",
           facecolor=star_color, edgecolors="white", linewidth=0.9,
           zorder=6, label="Test (evaluated at best epoch)")

# Handy labels
ax.set_xlabel("Epoch")
ax.set_ylabel("Loss (lower is better)")  # rename/invert if using a higher-is-better metric
ax.set_title("Chemprop (D-MPNN): Early stopping driven by validation", loc="left")

# Legend
leg = ax.legend(frameon=True, framealpha=0.96, facecolor="white", loc="upper right")
leg.get_frame().set_edgecolor("#e5e7eb")

fig.tight_layout()
# plt.savefig("chemprop_dmpnn_early_stopping_validation_driven.png", dpi=300, bbox_inches="tight")
# plt.savefig("chemprop_dmpnn_early_stopping_validation_driven.svg", bbox_inches="tight")
plt.show()