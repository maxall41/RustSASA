# $1 = Input
# $2 = Output
# $3 = Optional: --no-parallel to disable parallel processing

if [ "$3" = "--no-parallel" ]; then
    mkdir -p $2 && find $1 -type f | xargs -I {} bash -c 'uv run scripts/biopython_item.py "{}" "'$2'/$(basename "{}").json"'
else
    mkdir -p $2 && find $1 -type f | parallel uv run biopython_item.py {} $2/{/}.json
fi
