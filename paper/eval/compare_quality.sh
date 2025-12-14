sh scripts/get_data.sh
# run on e coli dataset
sh scripts/exec_freesasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_freesasa/
sh scripts/exec_rustsasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/
# run on freesasa dataset
sh scripts/exec_freesasa.sh ./data/pdb_files/ ./results/freesasa_dataset_proc=freesasa/
sh scripts/exec_rustsasa.sh ./data/pdb_files/ ./results/freesasa_dataset_proc=rustsasa/
# rename freesasa results to match file names (E. coli)
cd ./results/E_coli_freesasa
for f in *.pdb.json; do mv "$f" "${f/.pdb.json/.json}"; done
cd ../../
# rename freesasa results to match file names (freesasa dataset)
cd ./results/freesasa_dataset_proc=freesasa
for f in *.pdb.json; do mv "$f" "${f/.pdb.json/.json}"; done
cd ../../
# create figures
mkdir figures/
uv run scripts/quality_scatter.py ./results/E_coli_freesasa/ ./results/E_coli_RSASA/ "AlphaFold E. coli proteome" figures/sasa_chain_comparison_E_coli.pdf
uv run scripts/quality_scatter.py ./results/freesasa_dataset_proc=freesasa/ ./results/freesasa_dataset_proc=rustsasa/ "Freesasa dataset" figures/sasa_chain_comparison_freesasa_ds.pdf
uv run scripts/perf_graphs.py
