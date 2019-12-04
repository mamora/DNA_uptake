# DNA_uptake
Project data from Rosie Redfield laboratory. The goal of the project is to understand DNA uptake bias of Haemophilus influenzae. To fully understand the content of this repository you must read the paper that will be published about this project (DOI will be published as soon as it is published). 

Every folder in this repository contains an script, and information regarding the script, specific to a part of the analysis of the paper.

Raw sequencing files are linked to the Genbank project PRJNA387591 and will be available once the manuscript is published. Please see paper Table S2.  

You may request the authors for tables (dataframes) needed to run a particular script; but USS scores can be calculated with any genome sequence and uptake ratios could be calculated using the sequencing depth files linked to the bioproject number above. 

This folder are:

1. Sequencing_data_handling_scripts: This folder has shell scripts to calculate the samples sequencing depth from the raw sequencing files. For uptake samples, sequencing depth represents the level of uptake of donor genomic regions. A contamination correction takes place by adding a percentaje of random recipient read into the input donor file. Please refer to methods of the paper

2. Scoring the genome: This folder contain scripts that will "score" the NP and PittGG strains genome using a PWM (or PSSM) matrix derived from USS logos. This scores were used to create a list with genomic positions containing a USS-like sequence

3. Calculate_uptake_ratios: This folder scripts use the corrected sequencing depth from input and recovered periplasmic DNA to generate an uptake ratio using the normalize samples reads divided by the normalized input reads.

4. Predict uptake: Builds a deterministic model that predicts DNA uptake. Results from this model  will be used to analyze factors that affect DNA uptake.

5. Peak_symmetry: This script analyzed the symmetry of uptake peaks by calculating the mean peak shape of strong uptake peaks.

6. Make_sequence_motifs_of_uss_with_low_or_high_uptake: this folder has scripts that extract sequence of USS with low and high uptake and it build a motif logo out of those sequences

7. DNA_shape: This folder has a script that calculates the shape of DNA from a list of USS subset by level of uptake they promoted. It also calculate a Mann-Whitney-U statistical test to see differences in structural features between high and low uptake at each USS position. 


NOTES: To replicate this codes, you must change the paths within the scripts to fit where the files will be saved in your own computer. Also it is highly recomended that you work within a project folder in Rstudio or use the "here" package. 



