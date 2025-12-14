import json
import os
from collections import defaultdict

import fire
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr


def load_freesasa_data(file_path):
    """Load FreeSASA JSON data and extract chain areas."""
    with open(file_path) as f:
        data = json.load(f)

    chain_areas = {}
    for result in data["results"]:
        for structure in result["structure"]:
            for chain in structure["chains"]:
                chain_id = chain["label"]
                total_area = chain["area"]["total"]
                chain_areas[chain_id] = total_area

    return chain_areas


def load_rsasa_data(file_path):
    """Load RSASA JSON data and aggregate by chain."""
    with open(file_path) as f:
        data = json.load(f)

    chain_totals = defaultdict(float)
    for residue in data["Residue"]:
        chain_id = residue["chain_id"]
        value = residue["value"]
        chain_totals[chain_id] += value

    return dict(chain_totals)


def get_matching_files(dir1, dir2):
    """Get matching .json files from both directories."""
    files1 = set(f for f in os.listdir(dir1) if f.endswith(".json"))
    files2 = set(f for f in os.listdir(dir2) if f.endswith(".json"))
    return files1.intersection(files2)


def extract_matching_chain_values(freesasa_chains, rsasa_chains):
    """Extract values from matching chains."""
    freesasa_values = []
    rsasa_values = []

    for chain_id in freesasa_chains:
        if chain_id in rsasa_chains:
            freesasa_values.append(freesasa_chains[chain_id])
            rsasa_values.append(rsasa_chains[chain_id])

    return freesasa_values, rsasa_values


def main(freesasa_dir, rsasa_dir, plot_title, output_file):
    # Directory paths provided via CLI arguments

    # Check if directories exist
    if not os.path.exists(freesasa_dir):
        print(f"Directory {freesasa_dir} not found!")
        return
    if not os.path.exists(rsasa_dir):
        print(f"Directory {rsasa_dir} not found!")
        return

    # Get matching files
    matching_files = get_matching_files(freesasa_dir, rsasa_dir)

    if not matching_files:
        print("No matching .json files found in both directories!")
        return

    print(f"Found {len(matching_files)} matching files")

    # Collect all values from all files
    all_freesasa_values = []
    all_rsasa_values = []

    for filename in matching_files:
        freesasa_path = os.path.join(freesasa_dir, filename)
        rsasa_path = os.path.join(rsasa_dir, filename)

        try:
            freesasa_chains = load_freesasa_data(freesasa_path)
            rsasa_chains = load_rsasa_data(rsasa_path)

            values1, values2 = extract_matching_chain_values(
                freesasa_chains,
                rsasa_chains,
            )

            all_freesasa_values.extend(values1)
            all_rsasa_values.extend(values2)

            print(f"Processed {filename}: {len(values1)} matching chains")

        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue

    if not all_freesasa_values:
        print("No matching chains found!")
        return

    print(f"Total matching chains: {len(all_freesasa_values)}")

    # Convert to numpy arrays
    freesasa_values = np.array(all_freesasa_values)
    rsasa_values = np.array(all_rsasa_values)

    # Calculate Pearson correlation
    correlation, p_value = pearsonr(freesasa_values, rsasa_values)

    # Calculate additional statistics
    rmse = np.sqrt(np.mean((freesasa_values - rsasa_values) ** 2))
    mae = np.mean(np.abs(freesasa_values - rsasa_values))

    plt.figure(figsize=(10, 8))
    fontsize = 18
    linewidth = 2.5
    colors = ["#4477AA"]  # color-blind friendly blue from Paul Tol’s palette

    # Create axes
    ax = plt.subplot(111)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(linewidth)
    ax.spines["bottom"].set_linewidth(linewidth)
    ax.xaxis.set_tick_params(width=linewidth, length=8, direction="out")
    ax.yaxis.set_tick_params(width=linewidth, length=8, direction="out")

    ax.scatter(
        freesasa_values,
        rsasa_values,
        alpha=0.6,
        s=50,
        color=colors[0],
        edgecolor="none",
    )

    min_val = min(freesasa_values.min(), rsasa_values.min())
    max_val = max(freesasa_values.max(), rsasa_values.max())
    ax.plot([min_val, max_val], [min_val, max_val], color="gray", linestyle="--", linewidth=2)

    ax.set_xlabel("FreeSASA Chain Total", fontsize=fontsize + 2)
    ax.set_ylabel("RustSASA Chain Total", fontsize=fontsize + 2)
    ax.set_title(plot_title, fontsize=fontsize + 4, pad=15)

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
    ax.xaxis.get_offset_text().set_fontsize(fontsize - 2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize - 2)

    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_aspect("equal", adjustable="box")

    stats_text = f"n = {len(freesasa_values)}\nr = {correlation:.4f}\np = {p_value:.2e}\nRMSE = {rmse:.2f}"
    ax.text(
        0.05,
        0.95,
        stats_text,
        transform=ax.transAxes,
        va="top",
        fontsize=fontsize - 2,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight", dpi=300)

    print("\nResults:")
    print(f"Pearson correlation coefficient: {correlation:.4f}")
    print(f"P-value: {p_value:.2e}")
    print(f"Root Mean Square Error (RMSE): {rmse:.2f} Ų")
    print(f"Mean Absolute Error (MAE): {mae:.2f} Ų")
    print(f"Number of data points (chains): {len(freesasa_values)}")
    print("Plots saved as 'sasa_chain_comparison.png' and 'sasa_chain_comparison.pdf'")


if __name__ == "__main__":
    fire.Fire(main)
