#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
	cat <<EOF
Usage: $0 dir outdir
Assumes all files in dir are full of variants, and generates
matched random variants from each of these and spits them out in
outdir.

EOF
	exit 1
}

if [ $# -ne 2 ]; then
	usage
fi
in=$1
out=$2
max_jobs=50
split_count=50000
memory=25g
processors=1

outname=`basename $out`
logdir=~/sge_logs/silvafy/$outname
mkdir -pv $logdir
mkdir -pv $out/scripts


# Submit each file as a job, so long as it hasn't been completed
for file in $in/*; do
	if [ -d $file ]; then
		continue
	fi
	vcfname=`basename $file`
	if [[ -s "$out"/"$vcfname".random ]];  then
		echo "Already generated: $vcfname" >&2 
		continue
	fi
	mkdir -pv $out
	gen_script="$out/scripts/dispatch_gen_$vcfname.sh"
	cat > "$gen_script" <<EOF
#!/usr/bin/env bash
#$ -N "generate_$vcfname"
#$ -pe parallel 1
#$ -l h_vmem=10g
#$ -e $logdir
#$ -o $logdir
#$ -l hostname="supa*"

set -eu
set -o pipefail

source /dupa-filer/talf/silva-pipeline/silva_new/src/init.sh

temp=\$TMPDIR/"$vcfname".random
# Generate random mutations
/dupa-filer/talf/silva-pipeline/silva_new/src/src/input/synonymous.py generate $in/$vcfname -O /dupa-filer/talf/silva-pipeline/silva_new/src/data/refGene.pkl -G /dupa-filer/talf/silva-pipeline/silva_new/src/data/hg19.2bit --random --match-cpg > \$temp

# Two step move for safety
mv -v \$temp "$out"/"$vcfname".random.temp
mv -v "$out"/"$vcfname".random.temp "$out"/"$vcfname".random
EOF
qsub -S /bin/sh $gen_script
done
# Wait for them all to be done
slept="true"
while true; do
	slept="false"
	for file in $in/*; do
		if [ -d $file ]; then
			continue
		fi
		vcfname=`basename $file`
		if [[ ! -s "$out"/"$vcfname".random ]]; then
			echo $vcfname
			sleep 60 
			slept="true"
			break
		fi
	done
	if [ "$slept" = false ]; then
		echo here2
		break
	fi
	echo here
done

# Now submit silva jobs
for file in "$out"/*.random; do
	f=`basename $file`
	if [ -d $file ]; then
		continue
	fi
	outdir=$out/"$f"_processed
	mkdir -pv $outdir
	if [[ -s "$outdir"/0.scored ]]; then
		echo "VCF part already processed: $f" >&2
		continue
	fi

	script="$out/scripts/dispatch_$f.sh"
	cat > "$script" <<EOF
#!/usr/bin/env bash
#$ -N "$f"
#$ -pe parallel "$processors"
#$ -l h_vmem="$memory"
#$ -e $logdir
#$ -o $logdir
#$ -l hostname="supa*"

set -eu
set -o pipefail

# Our file is already filtered on AF, so we turn off AF filtering
export SILVA_AF_MAX=1
/dupa-filer/talf/silva-pipeline/silva_new/src/silva $outdir $file
EOF

	qsub -S /bin/sh $script
done
