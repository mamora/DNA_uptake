#align sequences to NP or GG genomes using BWA 
#          script                sample name       reference to use
bash alignerUptake3.sh pittgg250bp HiPittGG.fna
bash alignerUptake3.sh pittgg6000bp HiPittGG.fna


# mark duplicates and convert from sam files to bam files
#         script             input.name        reference      input.folder  output.folder
bash alignrest.sh pittgg250bp HiPittGG.fna sam bam
bash alignrest.sh pittgg250bp HiPittGG.fna sam bam