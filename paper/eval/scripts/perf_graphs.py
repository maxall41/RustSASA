import json
import os

import matplotlib.pyplot as plt

multi_data_fiiles = [
    "results/benches/bench_rustsasa.json",
    "results/benches/bench_freesasa.json",
    "results/benches/bench_biopython.json",
]
single_data_fiiles = [
    "results/benches/bench_rustsasa_single.json",
    "results/benches/bench_freesasa_single.json",
    "results/benches/bench_biopython_single.json",
]
singlethread_data_files = [
    "results/benches/bench_rustsasa_singlethread.json",
    "results/benches/bench_freesasa_singlethread_cpp.json",
    # "results/benches/bench_biopython_singlethread.json",
]


def create_panel(ax, json_files, title):
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
            if "rustsasa" in command or "rust-sasa" in command:
                label = "rust-sasa"
            elif "freesasa" in command or "sasa_batch" in command:
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

    # Determine if we should use milliseconds
    use_milliseconds = any(t < 1.0 for t in times)

    if use_milliseconds:
        times = [t * 1000 for t in times]
        errors = [e * 1000 for e in errors]
        time_unit = "milliseconds"
        time_suffix = "ms"
        decimal_places = 0
    else:
        time_unit = "seconds"
        time_suffix = "s"
        decimal_places = 1

    fontsize = 14
    linewidth = 2.0
    bar_colors = ["#33BBEE", "#EE7733", "#EE3377"]  # Paul Tol color-blind-safe palette

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

    ax.set_ylabel(f"Time ({time_unit})", fontsize=fontsize + 2)
    ax.set_title(title, fontsize=fontsize + 2, pad=10)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)

    for bar, err, val in zip(bars, errors, times):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + err + (0.02 * max(times)),
            f"{val:.{decimal_places}f}{time_suffix}",
            ha="center",
            va="bottom",
            fontsize=fontsize - 2,
        )

    max_height = max(times[i] + errors[i] for i in range(len(times)))
    ax.set_ylim(0, max_height * 1.15)
    ax.grid(True, axis="y", linestyle="--", alpha=0.3)


# Create a three-panel figure
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# Top left: multi_data_fiiles
ax1 = fig.add_subplot(gs[0, 0])
create_panel(
    ax1,
    multi_data_fiiles,
    r"$\mathbf{A.}$ Performance Comparison on E. coli Proteome",
)

# Top right: single_data_fiiles
ax2 = fig.add_subplot(gs[0, 1])
create_panel(
    ax2,
    single_data_fiiles,
    r"$\mathbf{B.}$ Performance Comparison on A0A385XJ53",
)

# Bottom: singlethread_data_files (spanning both columns)
ax3 = fig.add_subplot(gs[1, :])
create_panel(
    ax3,
    singlethread_data_files,
    r"$\mathbf{C.}$ Single-Threaded Performance on E. coli Proteome",
)

plt.savefig("figures/performance_comparison_combined.pdf", bbox_inches="tight", dpi=300)
