##################################################################################################
#         This script calculate the peaks shape of isolated USSs to assess symmetry              #
##################################################################################################

##############################################################
######   1. load samples and working directory        ########         
##############################################################

#  load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

###########################
#   load datasets         #
###########################
folder.name2 <- ("./datasets/final_datasets/") 

# uptake ratios
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")  

sigmodial_model<- fread(file = paste(folder.name2,
                                  "model_data/small_predicted_sigmodial_model_corrected_9.5.csv", 
                                  sep = ""))
# list of USS10
Up.USS.np.10.list<- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) 
# list of USS9.5
Up.USS.np.9.5.list<- fread(here::here("datasets/final_datasets","Up.USS.np.9.5.list.csv")) 

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")


##############################################################################
#           subset high uptake isolated USS and separate them by strand      #
##############################################################################

# calculate the distance to the next USS
next.uss.9.5<- sapply(Up.USS.np.9.5.list$keypos, FUN = dist.next.USS, 
                      USS.genome.c = Up.USS.np.9.5.list$keypos)
Up.USS.np.9.5.list$next.uss <- next.uss.9.5

# subset isolated USS
far<- which(Up.USS.np.9.5.list$next.uss >= 800)

isolated.uss.9.5.list<- Up.USS.np.9.5.list[far,]

# subset USS by strand **forward** strand = w and **reverse** strand = c

small.t<- isolated.uss.9.5.list[isolated.uss.9.5.list$strand == "w",]
small.t.c<- isolated.uss.9.5.list[isolated.uss.9.5.list$strand == "c",]

# subset only high scoring USS.
small_high_up<- small.t[which(small.t$keyup_small > 3), ] # forward

small_high_up.c<- small.t.c[which(small.t.c$keyup_small > 3), ] # reverse

##############################################################################
#           Build a matrix centered on USS +- 500 bp and calculate mean      #
#                               uptake peaks shape                           #
##############################################################################

#####################################
# for forward strand observed data  #
#####################################

# Next, I will built a matrix of each isolated USS with more than 3 uptake (rows) 
# and the uptake ratio of that USS minus/plus 50 bases (columns)
build_ratio_matrix = function( lrg=large , sml=small.c , radius=50 )
{
  nr = length (sml$USS.pos) 
  nc = 1+2*radius  
  out = matrix ( rep(NA,nr*nc) , ncol=nc , nrow=nr ) 
  for ( i in 1:nr ) 
  {
    out[ i ,] = lrg$ratio_short[max(1,sml$USS.pos[i]-radius):min(nrow(lrg),sml$USS.pos[i]+radius) ] 
  } 
  rownames(out) = as.character( sml$USS.pos ) 
  return( out ) 
}

#First concatenate the first 500 bases to the end to prevent errors and deal with circularity
uptake.f <- Uptake.ratio.np[1:500,]

# calculate distance from last USS in the genome to the position 1
dist.end <- (max(Up.USS.np.9.5.list$keypos) - length(Uptake.ratio.np$pos))  

# adjust distance to the end to circularize
uptake.c <- rbind(Uptake.ratio.np, uptake.f)  

# Subset the data that will be used to built a matrix of each uss +- 500 bases
large <- uptake.c[,c("pos","ratio_short")]
small <- small_high_up[,"USS.pos"]  

# Make vectors to build the matrix
al <- c(500:1,0,1:500) # make a vector of numbers from 500 to 0 to 500

# Make a vector of words "left, centre, right". 
# This words will be factors used to overlaplines in the following plot 
skew <- c(rep("Left",500),"centre",rep("Right",500)) 

# Calculate the average uptake ratio of each USS +- 500 bases for a given position in the USS. 
# I started with position 16 but this section be repeated for any positions in the motif.

# choose USS position to use as a centre. This number can be changed from 1 to 30
USS.pos <- 16 

# get all positions for each uss for the given USS.pos chosen
small.w <- small_high_up[,"USS.pos"]+(USS.pos -1) 

# make matrix centering on the USS.pos chosen
ratio_matrix.w = build_ratio_matrix(lrg=large, radius=500, sml = small.w) 

dim(ratio_matrix.w) # check dimensions. Should be 1001 columns (500+1+5000)

# calculate mean for each position
mean_col.w<- apply(ratio_matrix.w,2,mean, na.rm = TRUE) 

# calculate standard deviation for each position
sd_col.w<- apply(ratio_matrix.w,2,sd, na.rm = TRUE) 

flanks<- c(-500:500)

# Make a dataframe based on the matrix built before. GGplot won't work with a matrix
skew.test.w<- data.frame(flanks = flanks, pos = al, mean = mean_col.w, sd = sd_col.w,  align = skew) 

View(skew.test.w)

######################################
# for forward strand predicted data  #
######################################

model.w<- sigmodial_model[small.w$USS.pos]

# choose USS position to use as a centre. This number can be changed from 1 to 30
USS.pos<- 16 

small.m.w <- data.frame(USS.pos = model.w$pos)

# Subset the predicted data that will be used to built a matrix of each uss +- 500 bases
large<- data.frame(ratio_short = c(sigmodial_model$expected_uptake,sigmodial_model$expected_uptake[1:500]))

# make matrix centering on the USS.pos chosen
ratio_matrix.m.w = build_ratio_matrix(lrg=large, radius=500, sml = small.m.w) #make matrix centering on the USS.pos chosen

dim(ratio_matrix.m.w) # check dimensions. Should be 1001 columns (500+1+5000)

# calculate mean for each position
mean_col.w<- apply(ratio_matrix.m.w,2,mean, na.rm = TRUE) 

# make a dataframe based on the matrix built before. GGplot won't work with a matrix
skew.test.m.w<- data.frame(flanks = flanks, pos = al, mean = mean_col.w,  align = skew) 

###############################################################################################
#           Make plots to assess symnetry of isolated peaks based on matrix build before      #
###############################################################################################

#Important: Both lines should overlap perfectly when the highest point of uptake is found
p <- ggplot(aes(x = pos, y = mean, colour = align), data = skew.test.w) +
  geom_line() +
  scale_color_manual(values = c("Left" = "dark blue", "centre" = "black","Right" = "violetred3" ), name = "USS alignment") +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))+
  labs(x = "USS pos 15 +- 500 bases", y = "uptake ratio") +
  ggtitle(" Evaluate any skew of the curve of USS +- 500 bases in the forward strand") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = rel(1.2), face = "bold", vjust = 1.5))

p


p <- ggplot(aes(x = flanks, y = mean, colour = align), data = skew.test.w) +
  geom_line() +
  scale_color_manual(values = c("Left" = "dark blue", "centre" = "black","Right" = "violetred3" ), name = "USS alignment") +
  geom_line(aes(x = flanks, y = mean), colour = "dark green", alpha = 0.5, data = skew.test.m.w) +
  geom_errorbar(width=.1, aes(ymin=mean-sd, ymax=mean+sd),alpha = 1/10, colour="black") +
  scale_x_continuous(limits = c(-500,500), breaks = seq(-500, 500, by = 50), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0))+
  labs(x = "USS pos 15 +- 500 bases", y = "uptake ratio") +
  ggtitle(" Evaluate any skew of the curve of USS +- 500 bases in the forward strand") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = rel(1.2), face = "bold", vjust = 1.5))

p


p <- ggplot(aes(x = flanks, y = mean, colour = align), data = skew.test.w) +
  geom_line() +
  scale_color_manual(values = c("Left" = "dark blue", "centre" = "black","Right" = "violetred3" ), name = "USS alignment") +
  geom_line(aes(x = flanks, y = mean), colour = "dark green", alpha = 0.5, data = skew.test.m.w) +
  geom_errorbar(width=.1, aes(ymin=mean-sd, ymax=mean+sd),alpha = 1/10, colour="black") +
  scale_x_continuous(limits = c(-60,60), breaks = seq(-60, 60, by = 5), expand = c(0, 0))+
  scale_y_continuous(limits = c(2.5, 5.5), expand = c(0, 0))+
  labs(x = "USS pos 15 +- 500 bases", y = "uptake ratio") +
  ggtitle(" Evaluate any skew of the curve of USS +- 500 bases in the forward strand") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = rel(1.2), face = "bold", vjust = 1.5))

p

######################################
# for reverse strand observed data   #
######################################

large<- uptake.c[,c("pos","ratio_short")]

# choose USS position to use as a centre. This number can be changed from 1 to 30
USS.pos<- 16 

# get all positions for each uss for the given USS.pos chosen
small.c<- small_high_up.c[,"USS.pos"]+ 31 - USS.pos 

# make matrix centering on the USS.pos chosen
ratio_matrix.c = build_ratio_matrix(lrg=large, radius=500, sml = small.c) 

# check dimensions. Should be 1001 columns (500+1+5000)
dim(ratio_matrix.c) 

# calculate mean for each position
mean_col.c<- apply(ratio_matrix.c,2,mean, na.rm = TRUE) 

# calculate standard deviation for each position
sd_col.c<- apply(ratio_matrix.c,2,sd, na.rm = TRUE) 

# make a dataframe based on the matrix built before. GGplot won't work with a matrix
skew.test.c<- data.frame(flanks = flanks, pos = al, mean = mean_col.c, sd = sd_col.c,  align = skew) 

######################################
# for reverse strand predicted data  #
######################################

model.c<- sigmodial_model[small.c$USS.pos]

USS.pos<- 16 #chose USS position to use as a centre. This number can be changed from 1 to 30

small.m.c <- data.frame(USS.pos = model.w$pos)

large<- data.frame(ratio_short = c(sigmodial_model$expected_uptake,sigmodial_model$expected_uptake[1:500]))

# make matrix centering on the USS.pos chosen
ratio_matrix.m.c = build_ratio_matrix(lrg=large, radius=500, sml = small.m.c) 

dim(ratio_matrix.m.c) #check dimensions. Should be 1001 columns (500+1+5000)

mean_col.c<- apply(ratio_matrix.m.c,2,mean, na.rm = TRUE) #calculate mean for each position

# make a dataframe based on the matrix built before. GGplot won't work with a matrix
skew.test.m.c<- data.frame(flanks = flanks, pos = al, mean = mean_col.c,  align = skew) 

###############################################################################################
#           Make plots to assess symmetry of isolated peaks based on matrix build before      #
###############################################################################################


#Important: Both lines should overlap perfectly when the highest point of uptake is found
p <- ggplot(aes(x = pos, y = mean, colour = align), data = skew.test.c) +
  geom_line() +
  scale_color_manual(values = c("Left" = "dark blue", "centre" = "black","Right" = "violetred3" ), name = "USS alignment") +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))+
  labs(x = "USS pos 16 +- 500 bases", y = "uptake ratio") +
  ggtitle(" Evaluate any skew of the curve of USS +- 500 bases in the reverse strand") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = rel(1.2), face = "bold", vjust = 1.5))

p


p <- ggplot(aes(x = flanks, y = mean, colour = align), data = skew.test.c) +
  geom_line() +
  scale_color_manual(values = c("Left" = "dark blue", "centre" = "black","Right" = "violetred3" ), name = "USS alignment") +
  geom_line(aes(x = flanks, y = mean), colour = "dark green", alpha = 0.5, data = skew.test.m.w) +
  geom_errorbar(width=.1, aes(ymin=mean-sd, ymax=mean+sd),alpha = 1/10, colour="black") +
  scale_x_continuous(limits = c(-500,500), breaks = seq(-500, 500, by = 50), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0))+
  labs(x = "USS pos 16 +- 500 bases", y = "uptake ratio") +
  ggtitle(" Evaluate any skew of the curve of USS +- 500 bases in the reverse strand") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = rel(1.2), face = "bold", vjust = 1.5))

p


p <- ggplot(aes(x = flanks, y = mean, colour = align), data = skew.test.c) +
  geom_line() +
  scale_color_manual(values = c("Left" = "dark blue", "centre" = "black","Right" = "violetred3" ), name = "USS alignment") +
  geom_line(aes(x = flanks, y = mean), colour = "dark green", alpha = 0.5, data = skew.test.m.w) +
  geom_errorbar(width=.1, aes(ymin=mean-sd, ymax=mean+sd),alpha = 1/10, colour="black") +
  scale_x_continuous(limits = c(-60,60), breaks = seq(-60, 60, by = 5), expand = c(0, 0))+
  scale_y_continuous(limits = c(2.5, 5.5), expand = c(0, 0))+
  labs(x = "USS pos 16 +- 500 bases", y = "uptake ratio") +
  ggtitle(" Evaluate any skew of the curve of USS +- 500 bases in the reverse strand") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = rel(1.2), face = "bold", vjust = 1.5))

p