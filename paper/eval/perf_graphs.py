import matplotlib.pyplot as plt

# Data for plotting
libraries = ["rust-sasa", "freesasa", "biopython"]
times = [8.071, 54.914, 368.025]
errors = [0.361, 0.455, 51.156]

# --- Styling parameters ---
fontsize = 18
linewidth = 2.5
bar_colors = ["#33BBEE", "#EE7733", "#EE3377"]  # Paul Tol color-blind-safe palette

# --- Create figure and axes ---
fig, ax = plt.subplots(figsize=(10, 6))

# --- Axes styling ---
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_linewidth(linewidth)
ax.spines["bottom"].set_linewidth(linewidth)
ax.tick_params(axis="both", width=linewidth, length=8, direction="out")

# --- Bar chart (no outlines) ---
bars = ax.bar(
    libraries,
    times,
    yerr=errors,
    color=bar_colors,
    capsize=6,
)

# --- Labels ---
ax.set_ylabel("Time (seconds)", fontsize=fontsize + 2)
ax.set_title("Performance Comparison of SASA Libraries on E. coli Proteome", fontsize=fontsize + 4, pad=15)
ax.tick_params(axis="x", labelsize=fontsize)
ax.tick_params(axis="y", labelsize=fontsize)

# --- Optional: log scale (uncomment if large dynamic range) ---
# ax.set_yscale("log")

# --- Value labels above bars ---
for bar, err, val in zip(bars, errors, times):
    height = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        height + err + (0.02 * max(times)),
        f"{val:.1f}s",
        ha="center",
        va="bottom",
        fontsize=fontsize - 2,
    )

# --- Grid and limits ---
max_height = max(times[i] + errors[i] for i in range(len(times)))
ax.set_ylim(0, max_height * 1.15)
ax.grid(True, axis="y", linestyle="--", alpha=0.3)

# --- Layout and save ---
plt.tight_layout()
plt.savefig("figures/performance_comparison.pdf", bbox_inches="tight", dpi=300)
