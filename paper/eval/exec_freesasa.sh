# $1 = Input
# $2 = Output

mkdir -p $2 && find $1 -type f | parallel "freesasa --shrake-rupley -n 100 --format json {} > $2/{/}.json"
