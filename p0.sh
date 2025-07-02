#!/bin/bash

#Dependencies
#---
#cutdapt
#NGmerge

# set -x
set -e # Exit immediately if any of the following commands fail

READ1=$1
READ2=$2
OUTPATH=$3

# Sniff the thread count based on the number of logical CPUs, or the NUM_THREADS environment
# variable.
if [[ -z "${NUM_THREADS}" ]]; then
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	NT="$(nproc)"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
	NT="$(sysctl -n hw.logicalcpu)"
    else
        echo "Unknown OS ${OSTYPE}"
	exit 1
    fi
else
    NT=$NUM_THREADS
fi

if [[ "${NT}" -lt 2 ]]; then
    echo "Thread count ${NT} must be >= 2"
    exit 1
fi

mkdir -p $OUTPATH

#echo "Run FastQC"
#fastqc $READ1 -o $OUTPATH -t $NT --extract -q
#fastqc $READ2 -o $OUTPATH -t $NT --extract -q

echo "Trim illumina adapter sequences"
cutadapt -j $NT -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT -o $OUTPATH/R1_trimmed.fastq -p $OUTPATH/R2_trimmed.fastq $READ1 $READ2 --report minimal --json $OUTPATH/cutadapt.json > $OUTPATH/cutadapt.log

echo "Merge paired reads to one"
[[ -r $OUTPATH/merged.fastq ]] || NGmerge -1 $OUTPATH/R1_trimmed.fastq -2 $OUTPATH/R2_trimmed.fastq -o $OUTPATH/merged.fastq -n $NT -v -p 0.2 -m 10 2> $OUTPATH/merge.log

# Removed trimmed fastq files since no longer needed
rm $OUTPATH/R?_trimmed.fastq

{
    echo '{ "cutadapt": '
    cat $OUTPATH/cutadapt.json
    echo ','
    echo '"merge": {'
    grep Fragments $OUTPATH/merge.log | sed -e 's/.* /"input":/' -e 's/$/,/'
    grep Success $OUTPATH/merge.log | sed -e 's/.* /"output":/'
    echo '}'
    echo '}'
} > $OUTPATH/stats.json
