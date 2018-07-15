#!/bin/bash

prefix=$1 
refgen=$2 
inpdir=$3 
outdir=$4 

refpre=`basename $refgen .fna`

samblaster -i sam/$prefix.sam | sambamba view -h -f bam -S /dev/stdin | sambamba sort -o /dev/stdout /dev/stdin > $outdir/$prefix.$refpre.bam
sambamba index $outdir/$prefix.$refpre.bam
