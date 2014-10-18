#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
	cat <<EOF
Usage: $0 dir out

Combine the results of running silvafy on a vcf file.
You should pass the same dir you passed as outdir in silvafy.
EOF
	exit 1
}

if [ $# -ne 2 ]; then
	usage
fi

dir=$1
out=$2


export SILVA_AF_MAX=1

for subdir in $dir/*_processed; do
	/dupa-filer/talf/silva-pipeline/silva_new/src/silva-output.sh $subdir >> $2
done
