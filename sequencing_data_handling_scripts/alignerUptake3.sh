#!/bin/bash

prefix=$1 
refgen=$2 
inpdir=$3 
outdir=$4 
sam="@RG\tID:$prefix" 
refpre=`basename $refgen .fna`

/Linux/bin/bwa-0.7.15 mem -M -R $sam $refgen REMIX/$prefix.R1.fq.gz REMIX/$prefix.R2.fq.gz > sam/$prefix.sam


