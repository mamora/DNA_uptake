# make new concatenated genomes
cat ~/genomes/refs/rd.ncbi.fna ~/genomes/refs/np.ncbi.fna > ref/rdnp.fna
cat ~/genomes/refs/rd.ncbi.fna ~/genomes/refs/gg.ncbi.fna > ref/rdgg.fna

# Inidices
bwa index ref/rdnp.fna
bwa index ref/rdgg.fna

samtools faidx ref/rdnp.fna
samtools faidx ref/rdgg.fna

# simple beds
cat ref/rdnp.fna.fai | awk '{OFS="\t"; print $1, 0, $2}' > ref/rdnp.chr.bed
cat ref/rdgg.fna.fai | awk '{OFS="\t"; print $1, 0, $2}' > ref/rdgg.chr.bed