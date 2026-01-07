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


def process_dataset(freesasa_dir, rsasa_dir):
    """Process a single dataset and return the data arrays."""
    # Check if directories exist
    if not os.path.exists(freesasa_dir):
        print(f"Directory {freesasa_dir} not found!")
        return None, None
    if not os.path.exists(rsasa_dir):
        print(f"Directory {rsasa_dir} not found!")
        return None, None

    # Get matching files
    matching_files = get_matching_files(freesasa_dir, rsasa_dir)

    if not matching_files:
        print("No matching .json files found in both directories!")
        return None, None

    print(f"Found {len(matching_files)} matching files")

    # Collect all values from all files
    all_freesasa_values = []
    all_rsasa_values = []
    map_back = []

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
            map_back.append((len(all_freesasa_values), filename))

            print(f"Processed {filename}: {len(values1)} matching chains")

        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue

    if not all_freesasa_values:
        print("No matching chains found!")
        return None, None

    print(f"Total matching chains: {len(all_freesasa_values)}")

    dist = 0
    new_dist_id = None
    for i in range(len(all_freesasa_values)):
        new_dist = abs(all_rsasa_values[i] - all_freesasa_values[i])
        if new_dist > dist:
            dist = new_dist
            new_dist_id = i
    if new_dist_id is not None:
        print(map_back)
        print(f"Highest diff: {dist}: {new_dist_id}")

    # Convert to numpy arrays
    freesasa_values = np.array(all_freesasa_values)
    rsasa_values = np.array(all_rsasa_values)

    return freesasa_values, rsasa_values


def create_panel(ax, freesasa_values, rsasa_values, title):
    """Create a single scatter plot panel."""
    # Calculate Pearson correlation
    correlation, p_value = pearsonr(freesasa_values, rsasa_values)

    # Calculate additional statistics
    rmse = np.sqrt(np.mean((freesasa_values - rsasa_values) ** 2))

    fontsize = 14
    linewidth = 2.0
    colors = ["#4477AA"]  # color-blind friendly blue from Paul Tol's palette

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
    ax.set_title(title, fontsize=fontsize + 2, pad=10)

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

    return correlation, p_value, rmse


def main(
    freesasa_dir_ecoli,
    rsasa_dir_ecoli,
    freesasa_dir_freesasa,
    rsasa_dir_freesasa,
    output_file,
):
    # Process both datasets
    print("\nProcessing E. coli dataset...")
    ecoli_freesasa, ecoli_rsasa = process_dataset(freesasa_dir_ecoli, rsasa_dir_ecoli)

    print("\nProcessing Freesasa dataset...")
    freesasa_freesasa, freesasa_rsasa = process_dataset(freesasa_dir_freesasa, rsasa_dir_freesasa)

    if ecoli_freesasa is None or freesasa_freesasa is None:
        print("Error: Could not process one or both datasets")
        return

    # Create a two-panel figure
    fig = plt.figure(figsize=(16, 6))
    gs = fig.add_gridspec(1, 2, hspace=0.3, wspace=0.3)

    # Panel A: E. coli
    ax1 = fig.add_subplot(gs[0, 0])
    corr1, pval1, rmse1 = create_panel(
        ax1,
        ecoli_freesasa,
        ecoli_rsasa,
        r"$\mathbf{A.}$ AlphaFold E. coli proteome",
    )

    # Panel B: Freesasa dataset
    ax2 = fig.add_subplot(gs[0, 1])
    corr2, pval2, rmse2 = create_panel(ax2, freesasa_freesasa, freesasa_rsasa, r"$\mathbf{B.}$ FreeSASA dataset")

    plt.savefig(output_file, bbox_inches="tight", dpi=300)

    print("\n=== E. coli Results ===")
    print(f"Pearson correlation coefficient: {corr1:.4f}")
    print(f"P-value: {pval1:.2e}")
    print(f"Root Mean Square Error (RMSE): {rmse1:.2f} Ų")
    print(f"Number of data points (chains): {len(ecoli_freesasa)}")

    print("\n=== Freesasa Dataset Results ===")
    print(f"Pearson correlation coefficient: {corr2:.4f}")
    print(f"P-value: {pval2:.2e}")
    print(f"Root Mean Square Error (RMSE): {rmse2:.2f} Ų")
    print(f"Number of data points (chains): {len(freesasa_freesasa)}")

    print(f"\nCombined plot saved as '{output_file}'")


if __name__ == "__main__":
    fire.Fire(main)
