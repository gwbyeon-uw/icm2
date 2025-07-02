#!bin/bash

STRUCTFILE="structs.dot"
CSTART=16
CEND=75
EXPCL=3
NWORKERS=200
RDATFILE="Csde1_215_315.rdat"
OUTPREFIX="reeffit0"

reeffit --njobs 8 --structfile $STRUCTFILE --decompose motif --cstart $CSTART --cend $CEND --expcl $EXPCL --detailedplots --preparebootstrap --nworkers="$NWORKERS" --ntasks=1 $RDATFILE $OUTPREFIX

for file in `ls "$OUTPREFIX"bootstrap_worker*.sh`
do 
	cat bootstrap_slurm_header.txt $file > $file".run.sh"
	sbatch $file".run.sh"
done

