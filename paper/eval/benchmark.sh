 hyperfine --warmup 3 --runs 3 "sh exec_rustsasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/" --export-json ./results/bench_rustsasa.json
 hyperfine --warmup 3 --runs 3 "sh exec_freesasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/"  --export-json ./results/bench_freesasa.json
 hyperfine --warmup 3 --runs 3 "sh exec_biopython.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/"  --export-json ./results/bench_biopython.json
