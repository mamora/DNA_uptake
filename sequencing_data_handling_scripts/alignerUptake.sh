
prefix=$1 
refgen=$2 

refpre=`basename $refgen .fna` 


bwa mem -M -R '@RG\tID:$prefix' $refgen REMIX/$prefix.R1.fq.gz REMIX/$prefix.R2.fq.gz > re_bams/$prefix.sam
