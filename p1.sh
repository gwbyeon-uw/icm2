#!/bin/bash
# Dependencies: cutadapt, bowtie2, samtools, umicollapse, julia

#set -x
set -e # Exit immediately if any of the following commands fail

KEEPCLEAN=1   # Set to 1 to remove intermediate fastq files

INFASTQ=$1  # Input fastq file
OUTPATH=$2  # Folder in which to store output files
INDEX=$3    # Path to reference fasta file
SPLIT=$4    # Which split this is: 0-not split, 1-5' end, 2=3' end
SPLITSEQ=$5 # Seq to split on or blank to not split (always in same direction as sequence read p5->p7 )
UMILEN=$6   # Length of UMI (negative if on 3' end) (e.g. 14); 0 for no UMI
CONST5=$7   # Sequence of 5' constant region
CONST3=$8   # Sequence of 3' constant region
LTRIM=${#CONST5}
RTRIM=${#CONST3}

BINDIR="$(dirname $0)"

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

NEXT=$INFASTQ
if [[ $SPLIT == 1 ]]
then
    echo "Extracting split 1 concatenated reads"
    cutadapt -j $NT --discard-untrimmed -a $SPLITSEQ -o $OUTPATH/split.fastq $NEXT --report minimal -m 15 --json $OUTPATH/split.json > $OUTPATH/split.log
    NEXT=$OUTPATH/split.fastq
elif [[ $SPLIT == 2 ]]
then
    echo "Extracting split 2 concatenated reads"
    cutadapt -j $NT --discard-untrimmed -g $SPLITSEQ -o $OUTPATH/split.fastq $NEXT --report minimal -m 15 --json $OUTPATH/split.json > $OUTPATH/split.log
    NEXT=$OUTPATH/split.fastq
fi

if [[ $UMILEN != 0 ]]; then
    echo "UMI collapsing reads"
    umicollapse fastq -i $NEXT -o $OUTPATH/umicollapse.fastq > $OUTPATH/umicollapse.log
    [[ "$KEEPCLEAN" == 1 && "$NEXT" != "$INFASTQ" ]] && rm "$NEXT"
    NEXT=$OUTPATH/umicollapse.fastq
    echo "Remove UMI"
    cutadapt -j $NT -u $UMILEN -o $OUTPATH/umicut.fastq $NEXT -m 15 --json $OUTPATH/umicut.json > /dev/null # Cut off UMI
    [[ "$KEEPCLEAN" == 1 && "$NEXT" != "$INFASTQ" ]] && rm "$NEXT"
    NEXT=$OUTPATH/umicut.fastq
fi

# Trim constant regions
if [[ $CONST5 != "" || $CONST3 != "" ]]
then
    echo "Trimming constant regions"
    cutadapt -j $NT --discard-untrimmed --action=none --minimum-length 20 --rc -a '^'$CONST5...$CONST3'$' -o $OUTPATH/adaptercut.fastq $NEXT --report minimal --json $OUTPATH/primercut.json > $OUTPATH/primercut.log
    [[ "$KEEPCLEAN" == 1 && "$NEXT" != "$INFASTQ" ]] && rm "$NEXT"
    NEXT=$OUTPATH/adaptercut.fastq
fi

INDEX_NAME=$(basename $INDEX)
if [ -f "$OUTPATH/$INDEX_NAME"_idx"/$INDEX_NAME".1.bt2"" ]
then
    echo "Bowtie2 index exists: " $OUTPATH/$INDEX_NAME"_idx"/$INDEX_NAME
else
    echo "Building Bowtie2 index"
    # Build FASTA with constant regions from config added
    awk ' $0~/^>/ { header=$0; getline; print header; print "'$CONST5'"$0"'$CONST3'" } ' $INDEX > $OUTPATH/$INDEX_NAME
    mkdir -p $OUTPATH/$INDEX_NAME"_idx"
    bowtie2-build $OUTPATH/$INDEX_NAME $OUTPATH/$INDEX_NAME"_idx"/$INDEX_NAME > $OUTPATH/$INDEX_NAME"_idx.log"
fi

echo "Aligning reads"
bowtie2 -p $NT --end-to-end --very-sensitive -L 10 -N 1 --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 --no-unal -x $OUTPATH/$INDEX_NAME"_idx"/$INDEX_NAME -U $NEXT 2> $OUTPATH/bowtie.log | samtools view -bS - | samtools sort > $OUTPATH/out.bam
[[ "$KEEPCLEAN" == 1 && "$NEXT" != "$INFASTQ" ]] && rm "$NEXT"

samtools coverage $OUTPATH/out.bam > $OUTPATH/coverage.csv

echo "Making correlated mutation matrix" 
julia -t $NT --check-bounds=no "$BINDIR/m2matrix.jl" --bam $OUTPATH/out.bam --out $OUTPATH/m2matrix/ --plot true --minlength 0.9 --ltrim $LTRIM --rtrim $RTRIM --minmapq 42 > $OUTPATH/m2matrix.log
for d in $OUTPATH/m2matrix/*
do
    echo $d
    if [[ -r $d/align_1.dump ]]
    then
       echo Join $d/align_.dump
       cat $d/align_*.dump > $d/align.dump
       rm $d/align_*.dump
    fi
done

# Compile stats
{
echo "{"
if [[ $SPLIT != "0" ]]
then
    echo '"split":'
    cat $OUTPATH/split.json
    echo ','
fi

if [[ $UMILEN != 0 ]]; then
    echo '"umicollapse": {'
    grep 'input reads' $OUTPATH/umicollapse.log | sed -e 's/.*[ 	]/"inputreads":/' -e 's/$/,/'
    grep 'unique reads' $OUTPATH/umicollapse.log | sed -e 's/.*[ 	]/"uniquereads":/' -e 's/$/,/'
    grep 'deduplicating' $OUTPATH/umicollapse.log | sed -e 's/.*[ 	]/"deduppedreads":/'
    echo '},'
fi

echo '"primercut": '
cat $OUTPATH/primercut.json
echo ','

echo '"alignment": {'
     grep 'unpaired' $OUTPATH/bowtie.log | sed -e 's/^ *\([0-9]*\).*/"unpaired": \1,/'
     grep 'aligned 0' $OUTPATH/bowtie.log | sed -e 's/^ *\([0-9]*\).*/"unaligned": \1,/'
     grep 'aligned exactly 1' $OUTPATH/bowtie.log | sed -e 's/^ *\([0-9]*\).*/"aligned": \1,/'
     grep 'aligned >1' $OUTPATH/bowtie.log | sed -e 's/^ *\([0-9]*\).*/"ambiguous": \1/'
echo '},'

echo '"m2matrix": {'
#cat $OUTPATH/m2matrix.log
echo '}'
echo '}'
} > $OUTPATH/stats.json
