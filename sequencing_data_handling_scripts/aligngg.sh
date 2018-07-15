#!/bin/bash

# align create sam files
/Linux/bin/bwa-0.7.15 mem -M -R '@RG\tID:pittgg250bp' refs/HiPittGG.fna raw_reads/pittgg250bp.R1.fq.gz raw_reads/pittgg250bp.R2.fq.gz > sam/pittgg250.sam
/Linux/bin/bwa-0.7.15 mem -M -R '@RG\tID:pittgg6000bp' refs/HiPittGG.fna raw_reads/pittgg6000bp.R1.fq.gz raw_reads/pittgg250bp.R2.fq.gz > sam/pittgg6000.sam


# sam to bam
 samblaster -i sam/pittgg250bp.sam | sambamba view -h -f bam -S /dev/stdin | sambamba sort -o /dev/stdout /dev/stdin > bam/pittgg250.HiPittGG.bam
 
  samblaster -i sam/pittgg6000bp.sam | sambamba view -h -f bam -S /dev/stdin | sambamba sort -o /dev/stdout /dev/stdin > bam/pittgg6000.HiPittGG.bam

# index bam files  
  sambamba index bam/pittgg6000.HiPittGG.bam
  
  sambamba index bam/pittgg250.HiPittGG.bam
 
# make pilon corrected reference 
  java -Xmx16G -jar /Linux/bin/pilon.jar --genome refs/HiPittGG.fasta --frags bam/pittgg6000.HiPittGG.bam --output pilon_pittgg_pacBio

# make a depth per position table
  bedtools genomecov -d -ibam bam/pilon_pittgg250bp.pilon_pittgg_pacBio.bam > pilon_pittgg250bpbed
  
  bedtools genomecov -d -ibam bam/pilon_pittgg6000bp.pilon_pittgg_pacBio.bam > pilon_pittgg6000bp.bed
  
  
  