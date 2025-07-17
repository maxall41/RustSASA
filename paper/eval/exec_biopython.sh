# $1 = Input
# $2 = Output

mkdir -p $2 && find $1 -type f | parallel uv run biopython_item.py {} $2/{/}.json
