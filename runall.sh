#!/bin/bash
set -e
BINDIR=$(dirname $0)
for i in config/*.json
do
    BC=$(basename $i .json)
    bash $BINDIR/run.sh E06 ~/AAINGS6-fastq $i >> RUN-$BC.log 2>&1 &
done
