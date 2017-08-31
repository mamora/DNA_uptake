#####################################################################################
########       Scoring 86-028NP genome with Uptake-bias motif model       ###########
#####################################################################################


###########################################################
######   load samples and working directory        ########         
###########################################################


# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/" # Where these files are located

# Name the sites file, and genome reference file
fasta    <- "datasets/np.pilon.fasta" # Reference fasta to the NP genome

# Set working directory
setwd(whereami)

list.files(whereami) #see files in directory


# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("~/DNA_uptake/helper_functions/pssmFunctions1.R")

#read uptake pssm matrix, based on Mell et al 2012 uptake motif

uptake.pssm<- read.csv("~/DNA_uptake/datasets/degMotif_F.csv") # motif is loaded in a dataframe format, for the analysis a matrix format is needed.


View(uptake.pssm) #view the dataframe

str(uptake.pssm) #check the file

###########################################################
######   change format of the uptake motif  pssm   ########         
###########################################################


# load bases as valid characters 
four.chars <- chars[-length(chars)]  # gives chars takes all except last character "N" see chars in pssmFunctions1.R

char<- list(four.chars) #eliminate the n and make it a list

str(char) #check the list

uptake.pssm<- as.matrix(uptake.pssm[,2:31], dimnames = four.chars, byrow = TRUE) #make the dataframe into a matrix to use scoring function

str(uptake.pssm) #not ready dimension not right

dimnames(uptake.pssm) <- char #add right dimentions

str(uptake.pssm) #now dimensions are right

######################################################################################
######   disect scoring function to check if pssm was done properly           ########         
######################################################################################

##################### This is the core of the scoring matrix #######################
#Here I am disecting the function to generate two scoring matrixes used for scoring the forward pssnN and reverse strands pssnRC
sitelen <- dim(uptake.pssm)[2] # length of site
n       <- rep(0, dim(uptake.pssm)[2]) # add row for Ns
pssmN   <- rbind(uptake.pssm, n) # pssm with Ns for forward strand 
pssmRC  <- pssmN[c((length(chars)-1):1,length(chars)), #take row 4,3,2,1,5 for reverse strand
                 sitelen:1]  #take column 37:1


View(pssmRC) #order of rows should be t,g,c,a,n

View(pssmN) #order of row should be a,c,g,t,n


#everything fine? ok lets continue


###########################################################
######           load and score the genome         ########         
###########################################################


# Read in genome
genome <- read.genome(fasta) # read fasta with seqinr, then massage
lapply(genome, table)        # look okay? table(genome also works)


# score the genome (~6-7 minutes for ~2Mb)
### WARNING: SLOW ###
np.uptake.scores <- lapply(genome, scoreContig, scoremat=uptake.pssm, circle=T)


# For multi-contigs, usually should set circle=F
# Also requires changing subsequent code to handle lists (lapply)

# Only one contig?
np.uptake.scores <- np.uptake.scores[[1]]

# Check out the first few lines
head(np.uptake.scores, 50)

#Check out the score dataframe
str(np.uptake.scores)

#######################################
# Reduce to 1 genome length of scores #
#######################################

### Warning: Somewhat slow ###

# Get maximum score per pair of scores (w and c strands)
np.uptake.scores$max  <- apply(np.uptake.scores, 1, max)

# Get which strand had maximum score
np.uptake.scores$strand <- max.strand(np.uptake.scores)

View(np.uptake.scores) #check the score dataframe

write.csv(np.uptake.scores, file="datasets/np.uptake.scores.csv", quote=FALSE) #export uss scores to your working directory

np.uptake.scores<- read.csv("np.uptake.scores.csv") #read the file in case I am starting all over


np.uptake.scores[which.min(np.uptake.scores$w),] #check lowest score

np.uptake.scores[which.max(np.uptake.scores$w),] #check lowest score


np.uptake.scores[which.min(np.uptake.scores$c),] #check highest score

np.uptake.scores[which.max(np.uptake.scores$c),] #check highest score

#########################################################################################
#####################      Plot histograms of uptake scores         #####################
#########################################################################################

p <- ggplot() +
  geom_histogram(aes(x = w), binwidth = 0.2, colour = "black", data = np.uptake.scores) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.2e+5), expand = c(0, 0))+
  labs(x = "scores Watson strand") +
  ggtitle("Histogram of uptake scores for NP genome Watson strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","NP", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = w), binwidth = 0.2, colour = "black", data = np.uptake.scores) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 2e+3), expand = c(0, 0))+
  labs(x = "scores Watson strand") +
  ggtitle("Histogram of uptake scores for NP genome Watson strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","NP", "zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = c), binwidth = 0.2, colour = "black", data = np.uptake.scores) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.2e+5), expand = c(0, 0))+
  labs(x = "scores Crick strand") +
  ggtitle("Histogram of uptake scores for NP genome Crick strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","NP","c", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = c), binwidth = 0.2, colour = "black", data = np.uptake.scores) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 2e+3), expand = c(0, 0))+
  labs(x = "scores Crick strand") +
  ggtitle("Histogram of uptake scores for NP genome Crick strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","NP","c", "zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()

###########################################     histograms of both c and w together   ########################################


library(dplyr)
library(tidyr)

np.uptake.scores<- read.csv("./datasets/np.uptake.scores.csv")


scores.long.np<- np.uptake.scores[,c(1,2,3)] %>% tidyr::gather(strands, score, -X)

unique(scores.long.np$strands)

p <- ggplot() +
  geom_histogram(aes(x = score), binwidth = 0.2, colour = "black", data = scores.long.np) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1))+
  scale_y_continuous(limits = c(0, 3e+5), expand = c(0, 0))+
  labs(x = "USS scores") +
  ggtitle("Histogram of uptake scores for NP genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","NP","both", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = score), binwidth = 0.2, colour = "black", data = scores.long.np) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 3e+3), expand = c(0, 0))+
  labs(x = "USS scores") +
  ggtitle("Histogram of uptake scores for NP genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","NP","both","zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


### Same but for GG

gg.uptake.scores<- read.csv("./datasets/gg.uptake.scores.csv")


scores.long.gg<- gg.uptake.scores[,c(2,3,4)] %>% tidyr::gather(strands, score, -X)

unique(scores.long.gg$strands)

p <- ggplot() +
  geom_histogram(aes(x = score), binwidth = 0.2, colour = "black", data = scores.long.gg) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1))+
  scale_y_continuous(limits = c(0, 3e+5), expand = c(0, 0))+
  labs(x = "USS scores") +
  ggtitle("Histogram of uptake scores for PittGG genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","GG","both", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = score), binwidth = 0.2, colour = "black", data = scores.long.gg) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 3e+3), expand = c(0, 0))+
  labs(x = "USS scores") +
  ggtitle("Histogram of uptake scores for PittGG genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","GG","both","zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


################################        Normalize scores by genome size    ###########################

ph<- ggplot() +
  geom_histogram(aes(x = score, y=(..count.. * 1e+6)/1914386), binwidth = 0.2, colour = "black", data = scores.long.np)
  

p <- ggplot() +
  geom_histogram(aes(x = score, y=..count../1914386), binwidth = 0.2, colour = "black", data = scores.long.np) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1))+
  scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0))+
  labs(x = "USS scores", y  = "density") +
  ggtitle("Histogram of uptake scores for NP genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","NP","both","density","all", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = score, y=..count../1914386), binwidth = 0.2, colour = "black", data = scores.long.np) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 0.001), expand = c(0, 0))+
  labs(x = "USS scores") +
  ggtitle("Histogram of uptake scores for NP genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","NP","both","zoom","density", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()

ph1<- ggplot() +
  geom_histogram(aes(x = score, y=(..count.. * 1e+6)/1887050), binwidth = 0.2, colour = "black", data = scores.long.gg)


p1 <- ggplot() +
  geom_histogram(aes(x = score, y=..count../1887050), binwidth = 0.2, colour = "black", data = scores.long.gg) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1))+
  scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0))+
  labs(x = "USS scores") +
  ggtitle("Histogram of uptake scores for PittGG genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","GG","both","density","all", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p1)
dev.off()

p1 <- ggplot() +
  geom_histogram(aes(x = score, y=..count../1887050), binwidth = 0.2, colour = "black", data = scores.long.gg) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 0.001), expand = c(0, 0))+
  labs(x = "USS scores") +
  ggtitle("Histogram of uptake scores for PittGG genome both strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores","GG","both","zoom","density", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p1)
dev.off()


t<- ggplot_build(ph) 

count.data.np.m<- t$data[[1]]

tail(count.data.np.m)

t1<- ggplot_build(ph1) 

count.data.gg.m<- t1$data[[1]]


count.data.both<- data.frame(score = count.data.np$x, density.np = count.data.np$y, density.gg = count.data.gg$y )

count.data.both.m<- data.frame(score = count.data.np.m$x, count.m.np = count.data.np.m$y, count.m.gg = count.data.gg.m$y )


write.csv(count.data.both, "./datasets/count.data.both.csv")

count.data.both<- read.csv("./datasets/count.data.both.csv")

count.data.both$np_gg<- count.data.both$density.np/count.data.both$density.gg

count.data.both.m$np_gg<- count.data.both.m$count.m.np/count.data.both.m$count.m.gg


p2 <- ggplot() +
  geom_point(aes(x = score, y=np_gg), shape = 20, size = 3, colour = "black", data = count.data.both.m) +
  scale_x_continuous(limits = c(-2, 13), expand = c(0, 0))+
  scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0))+
  labs(x = "USS scores", y = "difference between counts.per.millon in NP/PittGG") +
  ggtitle("USS scores vs counts.per.millon of scores in NP/PittGG") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS scores/scores_counts.per.millon_np_gg", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p2)
dev.off()




#################################################################
#################Generate a USS list for NP    ##################
#################################################################

#Here we choose a cutoff of the minimum score that a sequence mst have to be consider a USS --> see manuscript methods

#cutoff = 10
Up.USS.np.10.w<- which(np.uptake.scores$w >= 10) #which positions belong to USS motif given a cuttoff score >= 10 (forward) 
Up.USS.np.10.c<- which(np.uptake.scores$c >= 10) #which positions belong to USS motif given a cuttoff score >= 10 (reverse)

#################################cutoff off 10######################################
w<- rep("w", times = length(Up.USS.np.10.w)) #create a vector of w elements

c<- rep("c", times = length(Up.USS.np.10.c)) #create a vector of c elements


w.scores<- np.uptake.scores$w[c(Up.USS.np.10.w)] #get the scores of forward strand

c.scores<- np.uptake.scores$c[c(Up.USS.np.10.c)] #get the scores of reverse strand


#make a dataframe list with uss positions, strand and score
Up.USS.np.10.list<- data.frame(strand = c(w,c), USS.pos = c(Up.USS.np.10.w,Up.USS.np.10.c), USS.score = c(w.scores, c.scores) )

str(Up.USS.np.10.list) 


#Reorder the list according to position
library(doBy)

Up.USS.np.10.list<- orderBy(~USS.pos+strand+USS.score, data = Up.USS.np.10.list)

write.csv(Up.USS.np.10.list, file="./datasets/Up.USS.np.10.list.csv", quote=FALSE) #save file





