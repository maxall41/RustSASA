# $1 = Input
# $2 = Output

mkdir $2
# ../../target/aarch64-apple-darwin/release/rust-sasa $1 $2 --format json
../../target/release/rust-sasa $1 $2 --format json
