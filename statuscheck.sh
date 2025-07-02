#!/bin/bash

# Usage: statuscheck E06  (folder containing Sxx)
for d in $@
do
    echo
    echo $d
    # Figure out stage of analysis
    stats="$d/stats.json"
    if [[ -r $d/cutadapt.json ]]
    then
       CUTIN=$(jq .read_counts.input $d/cutadapt.json)
       CUTOUT=$(jq .read_counts.output $d/cutadapt.json)
       echo " Cutadapt $CUTIN -> $CUTOUT $(($CUTOUT * 100 / $CUTIN))"
    fi
    if [[ -r $stats ]]
    then
       MERGEIN=$(jq .merge.input $stats)
       MERGEOUT=$(jq .merge.output $stats)
       FRAC=
       echo " Merge    $MERGEIN -> $MERGEOUT $(($MERGEOUT * 100 / (1+$MERGEIN)))"
    fi
    for s in $d/split*
    do
	[[ -d $s ]] || continue
	echo " $s"
	stats="$s/stats.json"
	if [[ -r $s/split.json ]]
	then
	    SPLITIN=$(jq .read_counts.input $s/split.json)
	    SPLITOUT=$(jq .read_counts.output $s/split.json)
	    echo "  Split   $SPLITIN -> $SPLITOUT $(($SPLITOUT * 100 / (1+$SPLITIN) ))"
	fi
	if [[ -r $stats ]]
	then
	    UMIIN=$(jq .umicollapse.inputreads $stats)
	    UMIOUT=$(jq .umicollapse.deduppedreads $stats)
	    echo "  UMI     $UMIIN -> $UMIOUT $(($UMIOUT * 100 / (1+$UMIIN)))"
	fi
	if [[ -r $s/primercut.json ]]
	then
	    CONSTIN=$(jq .read_counts.input $s/primercut.json)
	    CONSTOUT=$(jq .read_counts.output $s/primercut.json)
	    echo "  Const   $CONSTIN -> $CONSTOUT $(($CONSTOUT * 100 / (1+$CONSTIN)))"
	fi
	    
	if [[ -r $s/stats.json ]]
	then
	    BOWIN=$(jq .alignment.unpaired $stats)
	    BOWOUT=$(jq .alignment.aligned $stats)
	    echo "  Bowtie  $BOWIN -> $BOWOUT $(($BOWOUT * 100 / (1+$BOWIN))) $(($BOWOUT * 100 / $CUTIN))"
	fi
    done
done
