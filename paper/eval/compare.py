#!/usr/bin/env python3
"""Script to create scatter plot with Pearson correlation comparing
E_coli_freesasa and E_coli_RSASA SASA values.

FreeSASA files contain aggregate chain-level data, while RSASA files
contain individual residue data. This script aggregates RSASA data
by chain to compare with FreeSASA chain totals.
"""

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


def main(freesasa_dir, rsasa_dir, plot_title):
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

    # Create scatter plot
    plt.figure(figsize=(10, 8))
    plt.rcParams.update({"font.size": 20})
    plt.scatter(freesasa_values, rsasa_values, alpha=0.6, s=30)

    # Add 1:1 line for reference
    min_val = min(freesasa_values.min(), rsasa_values.min())
    max_val = max(freesasa_values.max(), rsasa_values.max())
    plt.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8, label="y=x")

    # Add labels and title
    plt.xlabel("FreeSASA Total")
    plt.ylabel("RustSASA Total")
    plt.title(
        plot_title,
    )
    plt.legend()

    # Add grid
    plt.grid(True, alpha=0.3)

    # Make plot square
    plt.axis("equal")

    # Add statistics text box
    stats_text = f"n = {len(freesasa_values)}\nr = {correlation:.4f}\np = {p_value:.2e}"
    plt.text(
        0.05,
        0.95,
        stats_text,
        transform=plt.gca().transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()

    # Save plot
    plt.savefig("sasa_chain_comparison.pdf", bbox_inches="tight")

    print("\nResults:")
    print(f"Pearson correlation coefficient: {correlation:.4f}")
    print(f"P-value: {p_value:.2e}")
    print(f"Root Mean Square Error (RMSE): {rmse:.2f} Ų")
    print(f"Mean Absolute Error (MAE): {mae:.2f} Ų")
    print(f"Number of data points (chains): {len(freesasa_values)}")
    print("Plots saved as 'sasa_chain_comparison.png' and 'sasa_chain_comparison.pdf'")


if __name__ == "__main__":
    fire.Fire(main)
