import json
import os

import matplotlib.pyplot as plt

# --- Load data from JSON files ---
json_files = [
    "results/bench_rustsasa.json",
    "results/bench_freesasa.json",
    "results/bench_biopython.json",
]

data_entries = []

for json_file in json_files:
    if not os.path.exists(json_file):
        print(f"Warning: {json_file} not found. Skipping.")
        continue

    with open(json_file) as f:
        data = json.load(f)
        result = data["results"][0]

        command = result["command"]
        mean_time = result["mean"]
        stddev = result["stddev"]

        label = "Unknown"
        if "rustsasa" in command:
            label = "rust-sasa"
        elif "freesasa" in command:
            label = "freesasa"
        elif "biopython" in command:
            label = "biopython"

        data_entries.append(
            {
                "label": label,
                "time": mean_time,
                "error": stddev,
            },
        )

# Sort order: rust-sasa, freesasa, biopython
order = {"rust-sasa": 0, "freesasa": 1, "biopython": 2}
data_entries.sort(key=lambda x: order.get(x["label"], 999))

libraries = [d["label"] for d in data_entries]
times = [d["time"] for d in data_entries]
errors = [d["error"] for d in data_entries]

fontsize = 18
linewidth = 2.5
bar_colors = ["#33BBEE", "#EE7733", "#EE3377"]  # Paul Tol color-blind-safe palette

fig, ax = plt.subplots(figsize=(10, 6))

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_linewidth(linewidth)
ax.spines["bottom"].set_linewidth(linewidth)
ax.tick_params(axis="both", width=linewidth, length=8, direction="out")

bars = ax.bar(
    libraries,
    times,
    yerr=errors,
    color=bar_colors,
    capsize=6,
)

ax.set_ylabel("Time (seconds)", fontsize=fontsize + 2)
ax.set_title("Performance Comparison of SASA Libraries on E. coli Proteome", fontsize=fontsize + 4, pad=15)
ax.tick_params(axis="x", labelsize=fontsize)
ax.tick_params(axis="y", labelsize=fontsize)

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

max_height = max(times[i] + errors[i] for i in range(len(times)))
ax.set_ylim(0, max_height * 1.15)
ax.grid(True, axis="y", linestyle="--", alpha=0.3)

plt.tight_layout()
plt.savefig("figures/performance_comparison.pdf", bbox_inches="tight", dpi=300)
