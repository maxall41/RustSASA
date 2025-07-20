# $1 = Input
# $2 = Output

mkdir -p $2 && find $1 -type f | parallel "freesasa --format json {} > $2/{/}.json"
