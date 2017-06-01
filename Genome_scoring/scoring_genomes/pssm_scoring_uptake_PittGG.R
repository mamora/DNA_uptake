#####################################################################################
########       Scoring PittGG genome with Uptake-bias motif model       ###########
#####################################################################################


###########################################################
######   load samples and working directory        ########         
###########################################################


# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/documentos phd/experiments/virtual experiments/Vex4/tutorialRjosh/Tutorial1/" # Where these files are located

# Name the sites file, and genome reference file
fasta    <- genome.PittGG # Reference fasta to the PittGG genome

# Set directory
setwd(whereami)

list.files(whereami) #see files in directory


# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("pssmFunctions1.R")

#read uptake pssm matrix, based on Mell et al 2012 uptake motif

uptake.pssm<- read.csv("datasets/degMotif_F.csv") # motif is loaded in a dataframe format, for the analysis a matrix format is needed.

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

# score the genome (~6-7 minutes for ~2Mb)
### WARNING: SLOW ###
PittGG.uptake.scores <- lapply(X = genome.PittGG, FUN = scoreContig, scoremat=uptake.pssm, circle=T)


# For multi-contigs, usually should set circle=F
# Also requires changing subsequent code to handle lists (lapply)

# Only one contig?
PittGG.uptake.scores <- PittGG.uptake.scores[[1]]

# Check out the first few lines
head(PittGG.uptake.scores, 50)

#Check out the score dataframe
str(PittGG.uptake.scores)

#######################################
# Reduce to 1 genome length of scores #
#######################################

### Warning: Somewhat slow ###

# Get maximum score per pair of scores (w and c strands)
PittGG.uptake.scores$max    <- apply(PittGG.uptake.scores, 1, max)

# Get which strand had maximum score
PittGG.uptake.scores$strand <- max.strand(PittGG.uptake.scores)


write.csv(PittGG.uptake.scores, file="C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Scoring_USS/Output/PittGG.uptake.scores.csv", quote=FALSE) #export uss scores



#################################################################
#################Generate a USS list for PittGG##################
#################################################################

#Here we choose a cutoff of the minimum score that a sequence mst have to be consider a USS --> see manuscript methods

#cutoff = 10
Up.USS.PittGG.10.w<- which(PittGG.uptake.scores$w >= 10) #which positions belong to USS motif given a cuttoff score >= 10 (forward) 
Up.USS.PittGG.10.c<- which(PittGG.uptake.scores$c >= 10) #which positions belong to USS motif given a cuttoff score >= 10 (reverse)

#################################cutoff off 10######################################
w<- rep("w", times = length(Up.USS.PittGG.10.w)) #create a vector of w elements

c<- rep("c", times = length(Up.USS.PittGG.10.c)) #create a vector of c elements


w.scores<- PittGG.uptake.scores$w[c(Up.USS.PittGG.10.w)] #get the scores of forward strand

c.scores<- PittGG.uptake.scores$c[c(Up.USS.PittGG.10.c)] #get the scores of reverse strand


#make a dataframe list with uss positions, strand and score
Up.USS.PittGG.10.list<- data.frame(strand = c(w,c), USS.pos = c(Up.USS.PittGG.10.w,Up.USS.PittGG.10.c), USS.score = c(w.scores, c.scores) )

str(Up.USS.PittGG.10.list) 


#Reorder the list according to position
library(doBy)

Up.USS.PittGG.10.list<- orderBy(~USS.pos+strand+USS.score, data = Up.USS.PittGG.10.list)

write.csv(Up.USS.PittGG.10.list, file="C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Scoring_USS/Output/Up.USS.PittGG.10.list.csv", quote=FALSE) #save file




