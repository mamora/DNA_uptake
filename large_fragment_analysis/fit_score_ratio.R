#######################################################################################
#                  USS score vs large fragment uptake ratio analysis                  #
#######################################################################################

# load packages
library(data.table)
library(here)
library(plyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(sicegar)

# load dataframes 
# load scores
np.USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")

# uptake ratios
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv") 

# load USS9.5 list
Up.USS.np.9.5.list<- fread(here::here("datasets/final_datasets","Up.USS.np.9.5.list.csv")) 

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")

# calculate distance to closest USS
next.uss.9.5<- sapply(Up.USS.np.9.5.list$keypos, FUN = dist.next.USS, 
                      USS.genome.c = Up.USS.np.9.5.list$keypos)

#####################################################################
#              define a window to do the analysis                   #
#####################################################################

#  To improve speed of functions and loops
n<- 1000
n1<- n/2
l<- length(Uptake.ratio.np$ratio_long)

#####################################################################
#              calculate mean uptake of 1kb window span             #
#####################################################################

# add bases to the end of the genome to circularize the genome
data.np<- c(Uptake.ratio.np$ratio_long[(l - (n1 - 1)):l], Uptake.ratio.np$ratio_long, Uptake.ratio.np$ratio_long[1:n1]) 

# calculate the mean uptake ratio of USS motif positions 
mean.ratio<- rollapply(data = data.np, width = (n + 1), align = "center", by = 1, FUN = mean) 

isTRUE(length(mean.ratio) == l)

Uptake.ratio.np$mean.ratio <- mean.ratio

mean.r<- Uptake.ratio.np$mean.ratio[Up.USS.np.9.5.list$USS.pos]

Up.USS.np.9.5.list$mean.ratio <- mean.r

#####################################################################
#                  fit a linear model and plot                      #
#####################################################################

# subset USS isolated by at least 1kb
test<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$next.uss >= 1000),]

fit <- lm(test$mean.ratio ~ test$USS.score, data = test)

summary(fit)

str(fit)

plot(test$USS.score, test$mean.ratio, axes = FALSE)
abline(fit, lwd=2)
axis(side = 1, at = c(9.5, 10, 10.5, 11, 11.5, 12, 12.5,13))
axis(side = 2, at = seq(0, 5, 0.2), las = 1) 
axis(side = 4, at = seq(0, 5, 0.2), las = 1) 


