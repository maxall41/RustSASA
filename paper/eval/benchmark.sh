 # hyperfine --warmup 3 --runs 3 "sh exec_rustsasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/" --export-json ./results/bench_rustsasa.json
 # hyperfine --warmup 3 --runs 3 "sh exec_freesasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/"  --export-json ./results/bench_freesasa.json
 # hyperfine --warmup 3 --runs 3 "sh exec_biopython.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/"  --export-json ./results/bench_biopython.json

hyperfine --warmup 3 --runs 25 -N "../../target/release/rust-sasa ./data/UP000000625_83333_ECOLI_v6/AF-A0A385XJ53-F1-model_v6.pdb ./results/test_out.pdb" --export-json ./results/bench_rustsasa_single.json
hyperfine --warmup 3 --runs 25 -N "freesasa --shrake-rupley -n 100 --format json ./data/UP000000625_83333_ECOLI_v6/AF-A0A385XJ53-F1-model_v6.pdb" --export-json ./results/bench_freesasa_single.json
hyperfine --warmup 3 --runs 25 -N "sh exec_biopython.sh ./data/UP000000625_83333_ECOLI_v6/AF-A0A385XJ53-F1-model_v6.pdb ./results/test_out" --export-json ./results/bench_biopython_single.json
