################################################
## Precalculate uptake initiation probability ##
################################################

# The goal of this script is to precalculate predicted uptake according to each USS score 
# according to the sigmoidal model of how uptake increase as a function of score and accroding to
# a linear model. 


# load USS scores 
USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")

library(data.table)

########################################################
#   calculate Uptake probability for sigmoidal model  ##
########################################################

# predict uptake for every score using sigmoidal model parameters
# sigmoidal equation based on the formula:   y ~ Imax/(1+expa(-a1*(x-tmid)) )
# where: “tmid” is score at half uptake, “Imax” is maximum uptake, 
# and “a1” is the slope at tmid.

# based on Caglar et al 2018
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5774301/


# Imax =  maximum uptake 3.77 (called "intensity" in Cablar paper)
Imax = 3.77
# a1 = slope at tmid 3.48
a1 = 3.48
# tmid =  score at half intensity  10.6 (called "time" in Caglar paper)
tmid = 10.6
# t = time (x parameter)

max_score<- round(max(USS.scores$max), digits = 2)
score <- seq(9, max_score, by = 0.01)
predicted_uptake<- c()
for (i in 1:length(score)){
  d0<- Imax/(1 + exp(-a1 *(score[i] - tmid)))
  predicted_uptake<- c(predicted_uptake, d0) 
}

sigmoidal_precalculated <- data.table(score, predicted_uptake)

write.csv(sigmoidal_precalculated, "sigmoidal_precalculated.csv")

########################################################
#   calculate Uptake probability for sigmoidal model  ##
########################################################

max_score<- round(max(USS.scores$max), digits = 3)
score <- seq(9, max_score, by = 0.01)

min_value<- round(min(score/max_score), digits = 3)
max_value<- round(max(score/max_score), digits = 3)

values<- round(score/max_score, digits = 3)

# predict uptake for every score using linear model where score of a USS is divided by the max score 
# using formula from the following reference:
# https://stats.stackexchange.com/questions/70801/how-to-normalize-data-to-0-1-range
linear_predicted_uptake <- round((values - min_value)/(max_value - min_value),digits = 3)

linear_precalculated <- data.table(score, linear_predicted_uptake)

write.csv(linear_precalculated, "linear_precalculated_start_at_zero.csv")


plot(score, linear_predicted_uptake)

