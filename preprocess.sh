#!/usr/bin/env bash

set -eu
set -o pipefail

source ${SILVA_PATH:-$(dirname $0)}/init.sh
export SILVA_AF_MAX=0.05
control=$SILVA_CONTROL.mat
src=$SILVA_PATH/src
data=$SILVA_PATH/data
modeldir=$SILVA_PATH/src/models/forest
traineddir=$SILVA_TRAINED

synonymous="$src/input/synonymous.py --genome=$data/hg19.fa.gz --genes=$data/refGene.ucsc.gz --cache-genes=$data/refGene.pkl"
gp1k="$src/input/1000gp.py $data/1000gp.refGene.vcf.gz $data/1000gp.refGene.pkl"

function usage {
    cat <<EOF
Usage: $0 OUTDIR VCF

Annotates rare, synonymous variants in VCF file.

Ranked results printed to stdout.

OUTDIR should be a directory specific to the VCF file.
Creates OUTDIR if it does not exit.

Existing files in OUTDIR are NOT overwritten, so this
script can be stopped and resumed.

Pretrained models found in: $traineddir
Will use TMPDIR='$TMPDIR' as needed.
EOF
    exit 1
}


if [[ $# -ne 2 ]]; then
    usage
fi
outdir="$(cd -P "$1"; pwd)"
vcf="$2"

if [[ ! -s "$vcf" ]]; then
    echo -e "Error: missing input file: $vcf\n" >&2
    usage
fi
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
# Filter to just synonymous variants
# Adds gene and transcript information

if [[ "$vcf" == *.pcoord ]]; then
    base=$(basename "$vcf" .pcoord)
    pcoord="--protein-coords"
else
    base=$(basename "$vcf" .vcf)
    pcoord=
fi
if [[ -n "$pcoord" ]]; then
    out=$base.flt
else
    out=$base.syn
fi
skip_if_exists $outdir/$out \
    && echo "Filtering for synonymous exonic variants..." >&2 \
    && $synonymous filter $pcoord "$vcf" | egrep -v "^Y\b" > $TMPDIR/$out \
    && echo "Removing variants on chromosome Y..." >&2 \
    && mv $TMPDIR/$out $outdir/$out \
    && echo "Left with $(grep -v '^#' $outdir/$out | wc -l) variants..." >&2
test -s $outdir/$out

# Possibly filter by 1000 Genomes Project AF
if [[ -z "$pcoord" && -n $SILVA_AF_MAX && -n $SILVA_AF_MIN ]]; then
    out=$base.af
    skip_if_exists $outdir/$out \
        && echo "Annotating with 1000 Genomes Project allele frequency..." >&2 \
        && cut -f 1,2,5 $outdir/$base.syn \
           | $gp1k \
           > $TMPDIR/$out \
        && mv $TMPDIR/$out $outdir/$out
    test -s $outdir/$out

    out=$base.flt
    skip_if_exists $outdir/$out \
        && echo "Filtering by allele frequency [$SILVA_AF_MIN, $SILVA_AF_MAX]..." >&2 \
        && paste $outdir/$base.af $outdir/$base.syn \
	   | awk -F"\t" '/^#/; /^[^#]/{OFS="\t"; if ($1 == ".") {$1 = 0} if ($1 >= '"${SILVA_AF_MIN}"' && $1 <= '"${SILVA_AF_MAX}"') {print $0}}' \
           | cut -f 2- \
           > $TMPDIR/$out \
        && mv $TMPDIR/$out $outdir/$out \
        && echo "Left with $(grep -v '^#' $outdir/$out | wc -l) variants..." >&2
    test -s $outdir/$out
elif [[ -z "$pcoord" ]]; then
    ln -s -T $outdir/$base.flt $outdir/$base.syn
fi
