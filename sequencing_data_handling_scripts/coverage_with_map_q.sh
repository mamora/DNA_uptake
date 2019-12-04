for bam in `ls up_bams/*I?.np.pilon.bam` 
do
  prefix=`basename $bam .bam`
  
  sambamba view -t 12 -f bam -F "mapping_quality > 0" $bam \
    | sambamba flagstat -t 12 /dev/stdin \
    > stats/$prefix.map.flag.txt

   bedtools coverage -d -abam $bam -b pilonrefs/np.pilon.chr.bed > depths/UP/$prefix.depth.bed
done


for bam in `ls up_bams/*I?.gg.pilon.bam`
do
  prefix=`basename $bam .bam`

  sambamba view -t 12 -f bam -F "mapping_quality > 0" $bam \
    | sambamba flagstat -t 12 /dev/stdin \
    > stats/$prefix.map.flag.txt
  bedtools coverage -d -abam $bam -b pilonrefs/gg.pilon.chr.bed > depths/UP/$prefix.depth.bed
done

