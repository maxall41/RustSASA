import json

import fire
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley


def calculate_sasa(pdb_file, output_file):
    # Parse PDB file
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    # Calculate SASA
    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    # Extract per-residue SASA values
    sasa_dict = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = f"{chain.id}_{residue.id[1]}"
                sasa_dict[res_id] = residue.sasa

    # Save to JSON
    with open(output_file, "w") as f:
        json.dump(sasa_dict, f, indent=2)


if __name__ == "__main__":
    fire.Fire(calculate_sasa)
