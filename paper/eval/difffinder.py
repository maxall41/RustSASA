#!/usr/bin/env python3
"""
Compare SASA values between FreeSASA and RustSASA output PDB files.
SASA values are stored in the B-factor column.
"""

from Bio.PDB import PDBParser
import numpy as np


def load_sasa_values(pdb_file):
    """Load SASA values (from B-factor column) for each residue."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    sasa_values = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                # Use first atom's B-factor as representative for the residue
                # (typically all atoms in a residue have the same SASA value)
                res_id = (chain.id, residue.id[1], residue.resname)
                for atom in residue:
                    sasa_values[res_id] = atom.bfactor
                    break  # Just take first atom

    return sasa_values


def compare_sasa(freesasa_file, rsasa_file, top_n=20):
    """Compare SASA values and find residues with largest differences."""
    print(f"Loading {freesasa_file}...")
    freesasa_values = load_sasa_values(freesasa_file)

    print(f"Loading {rsasa_file}...")
    rsasa_values = load_sasa_values(rsasa_file)

    # Calculate differences
    differences = []
    for res_id in freesasa_values:
        if res_id in rsasa_values:
            delta = abs(freesasa_values[res_id] - rsasa_values[res_id])
            differences.append({
                'chain': res_id[0],
                'resnum': res_id[1],
                'resname': res_id[2],
                'freesasa': freesasa_values[res_id],
                'rsasa': rsasa_values[res_id],
                'delta': delta
            })

    # Sort by delta (largest first)
    differences.sort(key=lambda x: x['delta'], reverse=True)

    # Print results
    print(f"\n{'='*80}")
    print(f"Top {top_n} residues with largest SASA differences:")
    print(f"{'='*80}")
    print(f"{'Rank':<6} {'Chain':<7} {'ResNum':<8} {'ResName':<8} {'FreeSASA':<12} {'RustSASA':<12} {'Delta':<10}")
    print(f"{'-'*80}")

    for i, diff in enumerate(differences[:top_n], 1):
        print(f"{i:<6} {diff['chain']:<7} {diff['resnum']:<8} {diff['resname']:<8} "
              f"{diff['freesasa']:<12.2f} {diff['rsasa']:<12.2f} {diff['delta']:<10.2f}")

    # Calculate statistics
    deltas = [d['delta'] for d in differences]
    print(f"\n{'='*80}")
    print(f"Statistics:")
    print(f"{'='*80}")
    print(f"Total residues compared: {len(differences)}")
    print(f"Mean absolute difference: {np.mean(deltas):.2f}")
    print(f"Median absolute difference: {np.median(deltas):.2f}")
    print(f"Max absolute difference: {np.max(deltas):.2f}")
    print(f"Min absolute difference: {np.min(deltas):.2f}")
    print(f"Std deviation: {np.std(deltas):.2f}")

    return differences


if __name__ == '__main__':
    freesasa_file = '1kqf_out_freesasa.pdb'
    rsasa_file = '1kqf_out_rsasa.pdb'

    differences = compare_sasa(freesasa_file, rsasa_file, top_n=20)
