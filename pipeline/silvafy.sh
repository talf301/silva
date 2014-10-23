#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
	cat <<EOF
Usage: $0 vcf outdir

Split a vcf file into chunks and process each with silva,
creating the resulting preprocessing subfolders. We assume
that this vcf has already been prefiltered for synonymous
and 1000 genomes variants, etc.
EOF
	exit 1
}

if [ $# -ne 2 ]; then
	usage
fi
vcf=$1
out=$2
max_jobs=50
split_count=50000
memory=25g
processors=1

function sge_wait {
	sleep_time=1  # seconds
	# Check SGE for number of jobs
	local sge_jobs="$(qstat | grep "EZR_" | wc -l)" || true
	while [[ $sge_jobs -ge $max_jobs ]]; do
		sleep $sleep_time
		sleep_time=$(expr $sleep_time "*" 2)
		sge_jobs="$(qstat | grep "EZR_" | wc -l)" || true
	done
}
outname=`basename $out`
logdir=~/sge_logs/silvafy/$outname
mkdir -pv $logdir
mkdir -pv $out/scripts

vcfname=`basename $vcf`

# Split the syn vcf into pieces
split -l $split_count "$vcf" $out/"$vcfname".new.

# Submit each file as a job, so long as it hasn't been completed
for file in $out/"$vcfname".new.*; do
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
source /dupa-filer/talf/silva-pipeline/silva_new/src/init.sh
export SILVA_AF_MAX=1
/dupa-filer/talf/silva-pipeline/silva_new/src/silva $outdir $file
EOF
	# Wait for space on the cluster
	sge_wait
	# Submit job
	qsub -S /bin/sh "$script"
done
