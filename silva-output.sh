#!/usr/bin/env bash

set -eu
set -o pipefail

source ${SILVA_PATH:-$(dirname $0)}/init.sh

control=$SILVA_CONTROL.mat
src=$SILVA_PATH/src
data=$SILVA_PATH/data
modeldir=$SILVA_PATH/src/models/forest
traineddir=$SILVA_TRAINED

synonymous="$src/input/synonymous.py --genome=$data/hg19.fa.gz --genes=$data/refGene.ucsc.gz --cache-genes=$data/refGene.pkl"
gp1k="$src/input/1000gp.py $data/1000gp.refGene.vcf.gz $data/1000gp.refGene.pkl"

function usage {
    cat <<EOF
Usage: $0 OUTDIR

This script only does the last part of silva, that is, given a directory
OUTDIR which has alredy had silva run, it prints the scored variants to
stdout.
EOF
    exit 1
}


if [[ $# -ne 1 ]]; then
    usage
fi
outdir="$(cd -P "$1"; pwd)"

if [[ ! -e $modeldir/test ]]; then
    echo "Could not find test script: $modeldir/test" >&2
    exit 1
fi

function skip_if_exists {
    local f="$1"
    if [[ -s "$f" ]]; then
	echo "Found existing: $f..." >&2
    fi
    test ! -s "$f"
}

init_message "$0" "$@"

mkdir -pv $outdir

### End of code originally in silva-preprocess ###
### Code originally in silva-run is below ###

fltfile="$outdir/$base.flt"
test -s "$fltfile"
	
matfile="$outdir/$base.input"
test -s "$matfile"

for modelfile in $traineddir/*.model; do
    pushd $modeldir > /dev/null
    if [[ ! -s $modelfile ]]; then
	echo "Error: could not find saved model: $modelfile" >&2
	exit 1
    fi
    
    # Run saved model on MAT file and created score file
    model=$(basename $modelfile .model)
    scorefile=$model.scored
    if [[ ! -s $outdir/$scorefile ]]; then
	echo "Running model $model..." >&2
	./test $modelfile $matfile | cut -f 1 > $outdir/.$scorefile \
	    && mv $outdir/.$scorefile $outdir/$scorefile
	test -e $outdir/$scorefile
    fi
    popd > /dev/null
done

# Print scored examples to stdout
echo -e "\nPrinting scored variants to stdout..." >&2
$SILVA_PATH/src/util/summarize_scores.py $fltfile $outdir/*.scored


echo "$0: SUCCESS" >&2
