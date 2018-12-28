#####################
###   Goal  #########
#####################
################################################################
# Fit a sigmodial equation to USS score vs uptake ratios       #
#      later include this in my prediction model               #
################################################################

# load packages
library(data.table)
library(here)
library(plyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(sicegar)

# load dataframes 
np.USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")

# uptake ratios
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv") #read uptake file 

# list of USSs, I made a list of sequences with score higher than 5 to assess 
# if very weak USS promote uptake  
Up.USS.np.5.list<- fread(here::here("datasets/final_datasets","Up.USS.np.5.list.csv")) #save file

# list of USS9.5
Up.USS.np.9.5.list<- fread(here::here("datasets/final_datasets","Up.USS.np.9.5.list.csv")) #save file

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")

# make a complete list of sequences with score > 5 < 9.5
Up.USS.np.5.list.e<- Up.USS.np.5.list[which(Up.USS.np.5.list$USS.score <  9.5),]

# calculate distance to the closest USS
next.uss.5<- sapply(Up.USS.np.5.list.e$keypos, FUN = dist.USS, 
                  USS.genome.c = Up.USS.np.9.5.list$keypos)
next.uss.9.5<- sapply(Up.USS.np.9.5.list$keypos, FUN = dist.next.USS, 
                    USS.genome.c = Up.USS.np.9.5.list$keypos)

# remove not needed columns
Up.USS.np.9.5.list <- Up.USS.np.9.5.list[,3:8]
Up.USS.np.5.list.e$V1 <- NULL

Up.USS.np.5.list.e$next.uss <- next.uss.5
Up.USS.np.9.5.list$next.uss <- next.uss.9.5

# subset isolated USS
far<- which(Up.USS.np.9.5.list$next.uss >= 800)
far.uss<- Up.USS.np.9.5.list[far,]
far2<- which(Up.USS.np.5.list.e$next.uss >= 800)
far.uss2<- Up.USS.np.5.list.e[far2,]

# join both USS lists
Up.USS.np<- rbind(far.uss2, far.uss)

#####################################################################
#        Calculate mean uptake ratio per all USS positions          #
#####################################################################

# circularize the genome
n<- 30
n1<- n/2
l<- length(Uptake.ratio.np$ratio_short)
# add bases to the end of the genome to circularize the genome
data.np<- c(Uptake.ratio.np$ratio_short[(l - (n1 - 1)):l], Uptake.ratio.np$ratio_short, Uptake.ratio.np$ratio_short[1:n1]) 


# calculate the mean uptake ratio of USS motif positions 
mean.ratio<- rollapply(data = data.np, width = (n + 1), align = "center", by = 1, FUN = mean) # calculate mean ratio of 31 positions left of each pos

isTRUE(length(mean.ratio) == l)

Uptake.ratio.np$mean.ratio <- mean.ratio

# order USS by position
Up.USS.np<- dplyr::arrange(Up.USS.np, USS.pos)

# remove first USSs that are not actually isolated since there is a USS close to end of genome--
Up.USS.np<- Up.USS.np[22:length(Up.USS.np$USS.pos),]

mean.r<- Uptake.ratio.np$mean.ratio[Up.USS.np$USS.pos]

Up.USS.np$mean.ratio <- mean.r

##########################################################
#        fit the sigmoidal function to the data          #
##########################################################

# sicegar package requieres variables names "time" and "intensity"
data.fit <- data.frame(time = Up.USS.np$USS.score, intensity = Up.USS.np$mean.ratio)

# check is data is in thr right format
dataCheck(data.fit)

# fit sigmoidal model
model.fit<- fitAndCategorize(data.fit)

str(model.fit)

norm.data<- model.fit$normalizedInput
test<- model.fit$sigmoidalModel

# plot sigmoidal curve with data
t<- figureModelCurves(norm.data, test, fittedXmin = 9.5) + 
  scale_x_continuous(limits = c(9.5,13), breaks = seq(9.5 , 13, 0.5), expand = c(0, 0))

t

# extract model variables to calculate sigmoidal curve "y-hat" 
# Imax =  maximum intensity 3.53
Imax = 3.53
# a1 = slope at tmid 3.62
a1 = 3.62
# tmid =  time at half intensity  10.6
tmid = 10.6
# t = time (x parameter)

# equation according to:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5774301/
time <- seq(9, 12.5, by = 0.01)
d1<- c()
for (i in 1:length(time)){
  d0<- Imax/(1 + exp(-a1 *(time[i] - tmid)))
  d1<- c(d1, d0) 
}

str(d1)

# plot sigmoidal curve again bot using ggplot
m.l<- data.table(time, d1)
p <- ggplot() +
  geom_point(aes(x = time, y = intensity), shape = 20, size = 4, data = data.fit) +
  geom_line(aes(x = time, y = d1), colour = "red", data = m.l) +
  scale_x_continuous(limits = c(9.5, 13), breaks = seq(9.5 , 13, 0.5)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0, 0)) +
  labs(x = "max uss score", y = "mean uptake ratio") +
  ggtitle("max uptake scores vs mean uptake ratio over 31 pos") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))


p


# analyze which USS with scores between 9.5 -10 promote uptake higher tha 2-fold background "~ 0.2" 
all <- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$next.uss >= 800 & Up.USS.np.9.5.list$USS.score < 10),] 

low<- all[which(all$keyup_small < 0.2),]

up<- all[which(all$keyup_small >= 0.2),]

