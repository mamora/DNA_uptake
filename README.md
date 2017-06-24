# DNA_uptake
Project data from Rosie Redfield laboratory. The goal of the project is to understand all factors involved in explaining variation in DNA uptake across the genome of Haemophilus influenzae

Every folder in this repository contains an script, and information regarding the script, specific to a part of the analysis of the paper

Most folders have a sequence, which means that you must run one in other to make the rest work. 

This folder are:

1. Genome _scoring: This folder contain scripts that will "score" the NP and PittGG strains genome using a PSSM matrix derived from USS logos. This scores were used to create a list with genomic positions containing a USS-like sequence

2. Contamination corrections: This folder has scripts to correct potential contamination from Rd recipient reads in NP donor uptake samples.  

3. calculate_uptake_ratios: This folder scripts use the corrected reads from previous folder to generate an uptake ratio using the normalize samples reads divided by the normalized input reads.

4. dist_closest_uss: calculates the distance from one genomic position to the closest USS in the list

5. list_of_isolated_uss: this scripts generate a list of uss that are far away from each other by a given distance 

6. skew_of_uss: This folder scripts seek to evaluate which positions correspond the highest DNA uptake in each USS peak. This is important since the USS in the forward and reverse strands need to be aligned.

7. align_uss_by_highest_uptake: Once the position with the highest uptake in each USS was identified. That position will be representing each USS in the list 

8. peak_finder: This script seek to generate a list of uptake peaks in the short fragment samples

9. syntenic_region_analysis: This scripts take a genomic region from a colinear alignment made by MAUVE and generate uptake maps from both NP and PittGG samples.
 

 
Other folders have scripts that do not belong to the main sequence analysis, but they run analysis that were needed for specific figures 

Some of this folders are:

1. Bioanalyzer_input_folder: analyzes the fragment distribution of the sheared DNA

2.  Figure_1_logos: This folder constain scripts that uses a modified seqlogo functions to remake the USS motifs from previous genomic and uptake bias analysis

 