#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
	cat <<EOF
Usage: $0 vcf outdir

Run the full analysis pipeline on a new file.
This involves generating matched random mutations,
running silva on both the original and matched, and then
producing some analysis results.
EOF
exit 1
}

function waitfor {
# $1 is the directory we're waiting on
	while true; do
		slept=false
		for $sub in $1/*_processed; do
			if [[ ! -s $sub/0.scored ]]; then
				sleep 60
				slept=true
				break
			fi
		done
		if [[ ! $slept ]]; then
			return
		fi
	done
}
if [[ $# -ne 2 ]]; then
	usage
fi

vcf=`pwd`/$1
outdir=`pwd`/$2
vcfname=`basename "$vcf" .vcf`

mkdir -pv "$outdir"

# Run our preprocessing so we get something we can generate from
/dupa-filer/talf/silva-pipeline/silva_new/src/preprocess.sh $outdir $vcf
flt=`basename "$outdir"/*.flt`

mkdir -pv "$outdir"/scripts
logdir="~/sge_logs/silva_generate/$2"
mkdir -pv "$logdir"

script="$outdir/scripts/dispatch_generate.sh"
cat > "$script" << EOF
#!/usr/bin/env bash
#$ -N "generate"
#$ -pe parallel 1
#$ -l h_vmem=5g
#$ -e $logdir
#$ -e $logdir
#$ -l hostname="supa*"

source /dupa-filer/talf/silva-pipeline/silva/init.sh

temp=\$TMPDIR/"$vcfname".random
# Generate random mutations
/dupa-filer/talf/silva-pipeline/silva_new/src/src/input/synonymous.py generate $outdir/$flt -O /dupa-filer/talf/silva-pipeline/silva_new/src/data/refGene.pkl -G /dupa-filer/talf/silva-pipeline/silva_new/src/data/hg19.2bit --random --match-cpg > $temp

# Two step move for safety
mv -v \$temp "$outdir"/"$vcfname".random.temp
mv -v "$outdir"/"$vcfname".random.temp "$outdir"/"$vcfname".random
EOF

# qsub mutation generation
if [[ ! -s "$outdir"/"$vcfname".random ]]; then
	qsub -S /bin/sh $script
fi

# SILVAFY
/dupa-filer/talf/silva-pipeline/silva/pipeline/silvafy.sh "$outdir"/"$flt" "$outdir"/original

# If mutations file doesn't exist, wait 1 minute and then check again
while [[ ! -s "$outdir"/"$vcfname".random ]]; do
	sleep 60
done

# SILVAFY MORE
/dupa-filer/talf/silva-pipeline/silva/pipeline/silvafy.sh "$outdir"/"$vcfname".random "$outdir"/random

# Wait for both of them to be done
waitfor "$outdir"/original
waitfor "$outdir"/random

# Combine outputs
if [[ -s "$outdir"/original_results.txt ]]; then
	read -p "File $outdir/original_results.txt already exists. Would you like to overwrite it?" -n 1 -r
	if [[ $REPLY =~ ^[Yy]$ ]]; then
		rm "$outdir"/original_results.txt
	else
		echo "Appending to existing file..."
	fi

/dupa-filer/talf/silva-pipeline/silva/pipeline/combine_pieces.sh "$outdir"/original "$outdir"/original_results.txt

# Combine outputs
if [[ -s "$outdir"/random_results.txt ]]; then
	read -p "File $outdir/random_results.txt already exists. Would you like to overwrite it?" -n 1 -r
	if [[ $REPLY =~ ^[Yy]$ ]]; then
		rm "$outdir"/random_results.txt
	else
		echo "Appending to existing file..."
	fi

/dupa-filer/talf/silva-pipeline/silva/pipeline/combine_pieces.sh "$outdir"/random "$outdir"/random_results.txt

# Run analysis
