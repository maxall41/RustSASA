uv run compare.py E_coli_freesasa/ E_coli_RSASA/ "Chain-Level SASA Comparison: E. coli dataset"
uv run compare.py freesasa_dataset_proc=freesasa/ freesasa_dataset_proc=rustsasa/ "Chain-Level SASA Comparison: Freesasa dataset"
uv run perf_graphs.py
