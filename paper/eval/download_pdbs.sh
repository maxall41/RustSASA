#!/bin/bash

# Input file containing PDB codes
INPUT_FILE="freesasa_dataset.txt"

# Create a directory for downloaded PDB files
mkdir -p pdb_files
cd pdb_files

# Read each PDB code from the file
while IFS= read -r pdb_code; do

  echo "Processing PDB code: $pdb_code"

  # Download the PDB file
  wget -q "https://files.rcsb.org/download/${pdb_code}.pdb"

  # Filter to keep only ATOM and HETATM lines
  grep -E "^(ATOM|HETATM)" "${pdb_code}.pdb" >"${pdb_code}_filtered.pdb"

  # Remove the original file and rename the filtered one
  rm "${pdb_code}.pdb"
  mv "${pdb_code}_filtered.pdb" "${pdb_code}.pdb"

done <"../$INPUT_FILE"
