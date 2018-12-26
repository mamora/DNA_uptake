# DNA_uptake
Project data from Rosie Redfield laboratory. The goal of the project is to understand all factors involved in explaining variation in DNA uptake across the genome of Haemophilus influenzae. To fully understand the content of this repository you must read the paper that will be published about this project around mid 2019 (hopefully). 

Every folder in this repository contains an script, and information regarding the script, specific to a part of the analysis of the paper.

Raw sequencing files are linked to a Genbank project nuber and will be available once the manuscript is published. Please see paper Table S1.  

Most folders are part of a sequence, which means that you must run one in other to make the rest work. Alternatively, you may request the authors for tables (dataframes) needed to run a particular script. 

This folder are:

1. sequencing_data_handling_scripts: This folder has shell scripts to calculate the samples sequencing depth from the raw sequencing files. For uptake samples, sequencing depth represents the level of uptake of donor genomic regions. A contamination correction takes place by adding a percentaje of random recipient read into the input donor file. Please refer to methods of the paper

2. scoring the genome: This folder contain scripts that will "score" the NP and PittGG strains genome using a PSSM matrix derived from USS logos. This scores were used to create a list with genomic positions containing a USS-like sequence

3. calculate_uptake_ratios: This folder scripts use the corrected reads from te first folder to generate an uptake ratio using the normalize samples reads divided by the normalized input reads.

4. dist_closest_uss: calculates the distance from one genomic position to the closest USS in the list and gets a list of USS "isolated" from ther USS by a given distance. In the paper we used 800bp

6. skew_of_uss: This folder scripts seek to evaluate which positions correspond the highest DNA uptake in each USS peak. This is important since the USS in the forward and reverse strands need to be aligned.

7. align_uss_by_highest_uptake: Once the position with the highest uptake in each USS was identified. That position will be representing each USS in the list 

9. syntenic_region_analysis: This scripts take a genomic region from a colinear alignment made by MAUVE and generate uptake maps from both NP and PittGG samples.
 

 Other folders have scripts that do not belong to the main sequence analysis, but they run analysis that were needed for specific figures 

Some of this folders are:

1. Bioanalyzer_input_folder: analyzes the fragment distribution of the sheared DNA

2.  Figure_1_logos: This folder constain scripts that uses a modified seqlogo functions to remake the USS motifs from previous genomic and uptake bias analysis

 