#!/bin/sh
# run.sh OUTROOT FASTQDIR [config ...]
# Config files listed on command line control processing of fastq files in FASTQDIR putting results in OUTROOT

set -e # Exit immediately if any of the following commands fail
if [[ ! -d $2 || ! -r $3 ]]
then
    echo Usage: $0 OUTROOT FASTQDIR [config ...]
    exit 1
fi

BINDIR=$(dirname $0)
OUTROOT=$1
shift
FASTQDIR=$1
shift

# Loop over the config files
for f in $*
do
    BC=$(jq -r '.barcode' $f)
    echo Processing barcode $BC
    # Destination for the analysis of this barcode
    OUTDIR=$OUTROOT/S${BC}
    # Check if the fastqs have already been merged
    if [ ! -r $OUTDIR/merged.fastq ]
    then
	# Process the fastq files to create merged.fastq (joined and adapters removed)
	$BINDIR/p0.sh $FASTQDIR/S${BC}_*R1*.fastq* $FASTQDIR/S${BC}_*R2*.fastq* $OUTDIR
    fi
    # Now use the merged.fastq as input to the next stage
    # For each split specified in config file
    nsplit=$(jq -r '.splits|length' $f)
    for i in $(seq 0 $(($nsplit-1)) )
    do
	SPLIT=$(jq -r ".splits[$i].splitnum" $f)
	MERGE=$(jq -r ".splits[$i].mergeseq" $f)
	UMILEN=$(jq -r ".splits[$i].umilength" $f)
	REF=$(jq -r ".splits[$i].referencefile" $f)
	MERGE=$(jq -r ".splits[$i].mergeseq" $f)
	CONST3=$(jq -r ".splits[$i].const3seq" $f)
	CONST5=$(jq -r ".splits[$i].const5seq" $f)
	       
	# New subdir for each split (i.e. for each different condition)
	base=$OUTDIR/split$SPLIT
	# Check if already processed
	if [ ! -r $base/out.bam ]
	then
	    # Locate the reference FASTA file
	    REF=$(dirname $f)/$(basename $REF)
	    if [[ ! -r $REF ]]
	    then
		echo Reference file, $REF, not found.
		exit 1
	    fi
	    # Run the rest of the pipeline (split,remove constant regions,align,compute m2seq matrix)
	    $BINDIR/p1.sh $OUTDIR/merged.fastq $base $REF $SPLIT $MERGE $UMILEN $CONST5 $CONST3
	fi
    done
done
