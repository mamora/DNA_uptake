#####################################################################################
########       Scoring NP genome with Uptake-bias motif model       ###########
#####################################################################################


###########################################################
######   load samples and working directory        ########         
###########################################################

# Name the sites file, and genome reference file
fasta    <- "./datasets/sequences/np.pilon.fasta" # Reference fasta to the NP genome

list.files(whereami) #see files in directory

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")

#read uptake pssm matrix, based on Mell et al 2012 uptake motif
# motif is loaded in a dataframe format, for the analysis a matrix format is needed.
uptake.pssm<- read.csv("./datasets/final_datasets/USS_scores/uptake_model/degMotif_F.csv") 
View(uptake.pssm) #view the dataframe
str(uptake.pssm) #check the file

###########################################################
######   change format of the uptake motif  pssm   ########         
###########################################################

# load bases as valid characters 
# gives chars takes all except last character "N" see chars in pssmFunctions1.R
four.chars <- chars[-length(chars)]  
char<- list(four.chars) #eliminate the n and make it a list
str(char) #check the list
#make the dataframe into a matrix to use scoring function
uptake.pssm<- as.matrix(uptake.pssm[,2:32], dimnames = four.chars, byrow = TRUE) 
str(uptake.pssm) #not ready dimension of the matrix are not right
dimnames(uptake.pssm) <- char #add the right dimentions
str(uptake.pssm) #now dimensions are right

###########################################################
######           load and score the genome         ########         
###########################################################

# Read in genome
genome <- read.genome(fasta) # read fasta with seqinr, then massage
lapply(genome, table)        # look okay? table(genome also works)



# score the genome (~6-7 minutes for ~2Mb)
### WARNING: SLOW ###
NP.uptake.scores <- lapply(X = genome, FUN = scoreContig, scoremat=uptake.pssm, circle=T)


# For multi-contigs, usually should set circle=F
# Also requires changing subsequent code to handle lists (lapply)

# Only one contig?
NP.uptake.scores <- NP.uptake.scores[[1]]

# Check out the first few lines
head(NP.uptake.scores, 50)

#Check out the score dataframe
str(NP.uptake.scores)

#######################################
# Reduce to 1 genome length of scores #
#######################################

### Warning: Somewhat slow ###

# Get maximum score per pair of scores (w and c strands)
NP.uptake.scores$max    <- apply(NP.uptake.scores, 1, max)

# Get which strand had maximum score
NP.uptake.scores$strand <- max.strand(NP.uptake.scores)

#export uss scores
write.csv(NP.uptake.scores, file="datasets/NP.USS.scores.pilon.corrected.csv", quote=FALSE) 

#########################################################################################
#####################      Plot histograms of uptake scores         #####################
#########################################################################################

p <- ggplot() +
  geom_histogram(aes(x = w), binwidth = 0.2, colour = "black", data = NP.uptake.scores) +
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
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores",
                  "NP", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = w), binwidth = 0.2, colour = "black", data = NP.uptake.scores) +
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
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores",
                  "NP", "zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = c), binwidth = 0.2, colour = "black", data = NP.uptake.scores) +
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
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores",
                  "NP","c", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = c), binwidth = 0.2, colour = "black", data = NP.uptake.scores) +
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
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores",
                  "NP","c", "zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()



#################################################################
#################Generate a USS list for NP   ##################
#################################################################

#Here we choose a cutoff of the minimum score that a sequence must have to be consider a USS

#cutoff = 9.5 # cutoff of 10 and 9.5 were used in the paper
#which positions belong to USS motif given a cuttoff score >= 9.5 (forward) 
Up.USS.NP.9.5.w<- which(NP.uptake.scores$w >= 9.5) 
#which positions belong to USS motif given a cuttoff score >= 9.5 (reverse)
Up.USS.NP.9.5.c<- which(NP.uptake.scores$c >= 9.5) 

#################################cutoff off 10######################################
w<- rep("w", times = length(Up.USS.NP.9.5.w)) #create a vector of w elements

c<- rep("c", times = length(Up.USS.NP.9.5.c)) #create a vector of c elements


w.scores<- NP.uptake.scores$w[c(Up.USS.NP.9.5.w)] #get the scores of forward strand

c.scores<- NP.uptake.scores$c[c(Up.USS.NP.9.5.c)] #get the scores of reverse strand


#make a dataframe list with uss positions, strand and score
Up.USS.NP.9.5.list<- data.frame(strand = c(w,c), 
                                    USS.pos = c(Up.USS.NP.9.5.w,Up.USS.NP.9.5.c), 
                                    USS.score = c(w.scores, c.scores) )

str(Up.USS.NP.9.5.list) 


#Reorder the list according to position
library(doBy)

Up.USS.NP.9.5.list<- orderBy(~USS.pos+strand+USS.score, data = Up.USS.NP.9.5.list)

write.csv(Up.USS.NP.9.5.list, 
          file="./datasets/final_datasets/USS_scores/uptake_model/Up.USS.NP.9.5.list.pilon.csv",
          quote=FALSE) #save file




