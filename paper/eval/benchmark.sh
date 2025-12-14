 hyperfine --warmup 3 --runs 3 "sh scripts/exec_rustsasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/" --export-json ./results/benches/bench_rustsasa.json
 hyperfine --warmup 3 --runs 3 "sh scripts/exec_freesasa.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/"  --export-json ./results/benches/bench_freesasa.json
 hyperfine --warmup 3 --runs 3 "sh scripts/exec_biopython.sh ./data/UP000000625_83333_ECOLI_v6/ ./results/E_coli_RSASA/"  --export-json ./results/benches/bench_biopython.json

hyperfine --warmup 3 --runs 25 -N "../../target/aarch64-apple-darwin/release/rust-sasa  ./data/UP000000625_83333_ECOLI_v6/AF-A0A385XJ53-F1-model_v6.pdb ./results/test_out.pdb" --export-json ./results/benches/bench_rustsasa_single.json
hyperfine --warmup 3 --runs 25 -N "freesasa --shrake-rupley -n 100 --format json ./data/UP000000625_83333_ECOLI_v6/AF-A0A385XJ53-F1-model_v6.pdb" --export-json ./results/benches/bench_freesasa_single.json
hyperfine --warmup 3 --runs 25 -N "sh scripts/exec_biopython.sh ./data/UP000000625_83333_ECOLI_v6/AF-A0A385XJ53-F1-model_v6.pdb ./results/test_out" --export-json ./results/benches/bench_biopython_single.json

hyperfine --warmup 3 --runs 3 "./sasa_batch ./data/UP000000625_83333_ECOLI_v6 ./out/" --export-json ./results/benches/bench_freesasa_singlethread_cpp.json
hyperfine --warmup 3 --runs 3 "../../target/aarch64-apple-darwin/release/rust-sasa  ./data/UP000000625_83333_ECOLI_v6/ ./out/ --format json --threads 1" --export-json ./results/benches/bench_rustsasa_singlethread.json
hyperfine --warmup 3 --runs 3 "sh scripts/exec_biopython.sh ./data/UP000000625_83333_ECOLI_v6/ ./out/ --no-parallel" --export-json ./results/benches/bench_biopython_singlethread.json
