#!/bin/bash
set -euo pipefail

mkdir -p data

# Get E coli AlphaFold dataset
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000625_83333_ECOLI_v6.tar
tar -xf UP000000625_83333_ECOLI_v6.tar -C data/
rm UP000000625_83333_ECOLI_v6.tar
gunzip data/UP000000625_83333_ECOLI_v6/*.gz
rm -rf data/UP000000625_83333_ECOLI_v6/*.cif

# Get Freesasa dataset PDB files

# Input file containing PDB codes
INPUT_FILE="./data/freesasa_dataset.txt"

# Create a directory for downloaded PDB files
mkdir -p data/pdb_files
cd data/pdb_files

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
