# DNA_uptake
Project data from Rosie Redfield laboratory. The goal of the project is to understand all factors involved in explaining variation in DNA uptake across the genome of Haemophilus influenzae. To fully understand the content of this repository you must read the paper that will be published about this project around mid 2019 (hopefully). 

Every folder in this repository contains an script, and information regarding the script, specific to a part of the analysis of the paper.

Raw sequencing files are linked to a Genbank project nuber and will be available once the manuscript is published. Please see paper Table S1.  

Most folders are part of a sequence, which means that you must run one in other to make the rest work. Alternatively, you may request the authors for tables (dataframes) needed to run a particular script. All of these datasets will be made available to public once the manuscript is ready to be submitted.   

This folder are:

1. Sequencing_data_handling_scripts: This folder has shell scripts to calculate the samples sequencing depth from the raw sequencing files. For uptake samples, sequencing depth represents the level of uptake of donor genomic regions. A contamination correction takes place by adding a percentaje of random recipient read into the input donor file. Please refer to methods of the paper

2. Scoring the genome: This folder contain scripts that will "score" the NP and PittGG strains genome using a PWM (or PSSM) matrix derived from USS logos. This scores were used to create a list with genomic positions containing a USS-like sequence

3. Calculate_uptake_ratios: This folder scripts use the corrected reads from te first folder to generate an uptake ratio using the normalize samples reads divided by the normalized input reads.

4. Build_fragment_class_matrix: This folder script takes fragment sizes as measured by bioanalyzer and subset fragment sizes into classes: 10bp/class for small fragment data and 200bp/class for large fragment data.

5. Ratio_vs_score: This folder scripts fit a sigmoidal function to understand hor uptake ratio increases as a function of USS score.

6. Predict uptake: Builds a deterministic model that predicts DNA uptake. Results from this model  will be used to analyze factors that affect DNA uptake.

7. Peak_symmetry: This script analyzed the symmetry of uptake peaks by calculating the mean peak shape of isolated strong peaks (> 800bp from the closest USS and with ratios > 3).

8. Motif_low_high_uptake: this fiolder has scripts that extract sequence of USS with low and high uptake and it build a motif logo out of those sequences

9. Plot_uptake maps: This folder has several scripts to plot uptake maps with standad deviation of replicates variation; as well as to plot a comparative uptake map of observed vs predicted uptake. This scripts were used to make most figures that include uptake maps.

10. Large_fragment_analysis: This folder has scripts to analyze uptake ratio variation in large fragment datasets. Basically, they assess the relevance of USS density and of USS score explaining uptake ratio variation.

11. DNA_shape_analysis: This folder has a script that calculates the shape of DNA from a list of USS subset by level of uptake they promoted. It also calculate a Mann-Whitney-U statistical test to see differences in structural features between high and low uptake at each USS position. 


NOTES: To replicate this codes, you must change the paths within the scripts to fit where the files will be saved in your own computer. Also it is highly recomended that you work within a project folder in Rstudio or use the "here" package. 

Scripts analysing PittGG dataset were basically identical to the ones posted here, just with a different dataset. For this reason I am not including them in the repository.



