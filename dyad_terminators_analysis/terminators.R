###########################################################################
#                       analysis of dyad terminators                      #
###########################################################################

# load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)

########################
#   load datasets      #
########################

# uptake ratios
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv") #read uptake file 

# list of USS10
Up.USS.np.10.list <- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) 

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")

# calculate distance to closest USS
next.uss<- sapply(Up.USS.np.10.list$keypos, FUN = dist.next.USS, 
                      USS.genome.c = Up.USS.np.10.list$keypos)

Up.USS.np.10.list$next.uss<- next.uss


##########################################################################################
#   build a matrix of pairs of USS close enough to each other to act as terminators      #
##########################################################################################

# subset only USS that are closer than 25bp
dye1<- which(Up.USS.np.10.list$next.uss <= 25)

# function to get even numbers
is.even <- function(x) x %% 2 == 0

# get USSs in pairs
index<- c(1:length(dye1))
even.index<- index[is.even(index) == TRUE]
odd.index<- index[is.even(index) == FALSE]

odd.dye<- dye1[odd.index]
even.dye<- dye1[even.index]

# make a dataset of pairs of USS with scores and uptake ratios
ter_pairs<- data.frame(USS1 = Up.USS.np.10.list$keypos[even.dye], USS2 = Up.USS.np.10.list$keypos[odd.dye],
    USS1_strand = Up.USS.np.10.list$strand[even.dye], USS2_strand = Up.USS.np.10.list$strand[odd.dye], 
    USS1_score = Up.USS.np.10.list$USS.score[even.dye] , USS2_score = Up.USS.np.10.list$USS.score[odd.dye])

# get the distance between pairs of USS
ter_pairs$separ<- abs(ter_pairs$USS2 - ter_pairs$USS1)

# calculate mean uptake ratio 15 bp upstream of pos 16 from the first USS to 15 bp downstream of
# pos 16 from the second USS
ratio<- c()    
for(i in 1:length(ter_pairs$USS1)){
ra<- mean(Uptake.ratio.np$ratio_short[(ter_pairs$USS1[i] - 15)]:Uptake.ratio.np$ratio_short[(ter_pairs$USS2[i] + 15)]) 
ratio<- c(ratio, ra)
}

ter_pairs$ratio<- ratio

# save in file
write.csv(ter_pairs, file = "./Marcelo_paper_figures_new_order/supplementary/suppl_4/terminators/terminators.csv")

# subset pairs of USS by orientation
w.c<- which(ter_pairs$USS1_strand == "w" & ter_pairs$USS2_strand == "c")

c.w<- which(ter_pairs$USS1_strand == "c" & ter_pairs$USS2_strand == "w")

w.w<- which(ter_pairs$USS1_strand == "w" & ter_pairs$USS2_strand == "w")

c.c<- which(ter_pairs$USS1_strand == "c" & ter_pairs$USS2_strand == "c")

# get the average USS score
average<- round((ter_pairs$USS1_score + ter_pairs$USS2_score)/2, digits = 1)

#####################################################################################
# calculate the sigmoidal curve from a USS score vs uptake ratio analysis. Figure 3 #
#####################################################################################


# Imax =  maximum intensity 3.53
Imax = 3.53
# a1 = slope at tmid 3.62
a1 = 3.62
# tmid =  time at half intensity  10.6
tmid = 10.6
# t = time (x parameter)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5774301/
time <- seq(9, 13, by = 0.01)
d1<- c()
for (i in 1:length(time)){
  d0<- Imax/(1 + exp(-a1 *(time[i] - tmid)))
  d1<- c(d1, d0) 
}


# get "y-hat" for the sigmoidal curve or the predicted uptake ratio  if the score/ratio relationship 
# follows the sigmoidal curve  

exp.up<- c()
for(i in 1:length(average)){
exp_averg_upt<- d1[which(time %in% abs(average[i]))] 
exp.up<- c(exp.up, exp_averg_upt)
}

ter_pairs$exp_upt_mean<- exp.up
ter_pairs$diff_exp_mean<- round(ter_pairs$ratio - ter_pairs$exp_upt_mean, digits = 2)

###########################################################################
# calculate which member of the USS pair has the highest and lowest score #
###########################################################################
# which member pair has the highest score
max<- c()
for(i in 1:length(ter_pairs$USS1)){
m<- max(ter_pairs$USS1_score[i], ter_pairs$USS2_score[i])
max<- c(max,m)
}

max<- round(max, digits = 1)

# calculate predicted uptake for a USS with the max pair USS
exp.up<- c()
for(i in 1:length(average)){
  exp_averg_upt<- d1[which(time %in% abs(max[i]))] 
  exp.up<- c(exp.up, exp_averg_upt)
}

ter_pairs$exp_upt_max <- exp.up
ter_pairs$diff_exp_max<- round(ter_pairs$ratio - ter_pairs$exp_upt_max, digits = 2)

# which member pair has the lowest score
min<- c()
for(i in 1:length(ter_pairs$USS1)){
  m<- min(ter_pairs$USS1_score[i], ter_pairs$USS2_score[i])
  min<- c(min,m)
}

min<- round(min, digits = 1)

# calculate predicted uptake for a USS with the lowest uptake in the pair-USS
exp.up<- c()
for(i in 1:length(average)){
  exp_averg_upt<- d1[which(time %in% abs(min[i]))] 
  exp.up<- c(exp.up, exp_averg_upt)
}

ter_pairs$exp_upt_min <- exp.up
ter_pairs$diff_exp_min<- round(ter_pairs$ratio - ter_pairs$exp_upt_min, digits = 2)

# calculate the difference between the observed uptake and the predicted by the sigmoidal model
l<- length(ter_pairs$diff_exp_min)
diff<- c(ter_pairs$diff_exp_mean, ter_pairs$diff_exp_min, ter_pairs$diff_exp_max)
samples<- c(rep("average diff", l), rep("min diff", l), rep("max diff", l))

# get the max and min scores of the USS pair
max<- c()
min<-c()
for(i in 1:length(ter_pairs$USS1)){
  mx<- max(ter_pairs$USS1_score[i],ter_pairs$USS2_score[i] )
  mn<- min(ter_pairs$USS1_score[i],ter_pairs$USS2_score[i] )
  max<- c(max, mx)
  min<- c(min, mn)
}

m.l<- data.frame(time, d1)
ter_pairs$max.score<- max
ter_pairs$min.score<- min
  

# add USS orientation to the dataframe
ter_pairs$orientation<- rep("NA", 124)
ter_pairs$orientation[w.c]<- "w.c"
ter_pairs$orientation[c.w]<- "c.w"

###############################################################################
#   plot USS terminators pairs USS score vs ratio subset by orientation       #
###############################################################################

p <- ggplot() +
  geom_point(aes(x = min.score, y = ratio, colour = orientation), shape = 20, size = 2, data = ter_pairs) +
  geom_line(aes(x = time, y = d1), colour = "red", data = m.l) +
  scale_x_continuous(limits = c(9.5, 13), breaks = seq(9.5 , 13, 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0, 0)) +
  labs(x = "min uss score of a pair of USS", y = "mean uptake ratio") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))


p

p <- ggplot() +
  geom_point(aes(x = max.score, y = ratio, colour = orientation), shape = 20, size = 2, data = ter_pairs) +
  geom_line(aes(x = time, y = d1), colour = "red", data = m.l) +
  scale_x_continuous(limits = c(9.5, 13), breaks = seq(9.5 , 13, 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0, 0)) +
  labs(x = "max uss score of a pair of USS", y = "mean uptake ratio") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))


p


##############################################################################################
#   plot USS terminators pairs USS score vs ratio subset by distance betwen each other       #
##############################################################################################

# subset piar by dustance between each other
ter_pairs$fact<- rep("NA", 124)
ter_pairs$fact[which(ter_pairs$separ == 0)]<- "zero"
ter_pairs$fact[which(ter_pairs$separ == 1)]<- "one"
ter_pairs$fact[which(ter_pairs$separ == 2)]<- "two"
ter_pairs$fact[which(ter_pairs$separ == 3)]<- "three"
ter_pairs$fact[which(ter_pairs$separ >= 9)]<- "more than 9"

unique(ter_pairs$separ)

p <- ggplot() +
  geom_point(aes(x = min.score, y = ratio, colour = fact), shape = 20, size = 2, data = ter_pairs) +
  geom_line(aes(x = time, y = d1), colour = "red", data = m.l) +
  scale_color_manual(values = c("zero" = "green1","one" = "gold2", "two" = "darkorange1", "three" = "pink", "more than 9" = "red3"), name = "spacing between USSs") +
  scale_x_continuous(limits = c(9.5, 13), breaks = seq(9.5 , 13, 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0, 0)) +
  labs(x = "min uss score of a pair of USS", y = "mean uptake ratio") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))


p

p <- ggplot() +
  geom_point(aes(x = max.score, y = ratio, colour = fact), shape = 20, size = 2, data = ter_pairs) +
  geom_line(aes(x = time, y = d1), colour = "red", data = m.l) +
  scale_color_manual(values = c("zero" = "green1","one" = "gold2", "two" = "darkorange1", "three" = "pink", "more than 9" = "red3"), name = "spacing between USSs") +
  scale_x_continuous(limits = c(9.5, 13), breaks = seq(9.5 , 13, 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0, 0)) +
  labs(x = "max uss score of a pair of USS", y = "mean uptake ratio") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))


p

one<- ter_pairs[which(ter_pairs$separ >= 9),]
