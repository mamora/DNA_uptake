#####################################################################################
########       Scoring 86-028NP genome with Genomic motif model       ###########
#####################################################################################

###########################################################
######   load samples and working directory        ########         
###########################################################


# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/documentos phd/experiments/virtual experiments/Vex4/tutorialRjosh/Tutorial1/" # Where these files are located

# Name the sites file, and genome reference file
fasta    <- "np.fa" # Reference fasta to the NP genome
sites    <- "Hin2206sites.txt" # From Rosie's Gibbs sampling. check: Maughan et al. [2010] 

# Set working directory
setwd(whereami)

list.files(whereami) #see files in directory


# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("pssmFunctions1.R")


##################################################################
################## PART 1: Generate Scoring Matrix ###############
##################################################################

#############
# Make PSSM #
#############

sitemat <- sites2matrix(sites)    # set of sites
pcm     <- matrix2pcm(sitemat)    # count matrix
pfm     <- pcm2pfm(pcm, pseudo=0) # frequency matrix
pssm    <- pfm2pssm(pfm)          # scoring matrix


pssm <-  pssm[1:4,3:32] #remove the last positions that are not important
pfm <-   pfm[1:4,3:32]

ic      <- colSums(pssm)          # information content per position
wblg    <- pfm2wblg(pfm)          # matrix like weblogo output

# Determine the consensus sequence
consensus.num  <- apply(pfm, MARGIN=2, which.max) #  highest frequency base per position
consensus.site <- sapply(consensus.num, function(a) toupper(chars[a])) # convert to a letter


###################
# Dissect Weblogo #
###################

# Base barplots of pfm, ic, and wblg
wblgcols<- c("green", "blue", "gold", "red") # emulate weblogo base colors for ACGT

# Make a plot

barplot(wblg, names.arg=consensus.site, col=wblgcols, main="Genomic USS logo ", ylab="bits", ylim = c(0,2))
  axis(2, at = c(0,0.5,1,2))



#save file
write.csv(pssm, file="pssm.fixed.csv", quote=FALSE) #save file

write.csv(pfm, file="pfm.fixed.csv", quote=FALSE) #save file


##################################################################
############### PART 2: Score a genome with matrix ###############
##################################################################

################################
# Score genome on both strands #
################################

# Read in genome
genome <- read.genome(fasta) # read fasta with seqinr, then massage
lapply(genome, table)        # look okay? table(genome also works)


# score the genome (~6-7 minutes for ~2Mb)
### WARNING: SLOW ###
np.uss.scores <- lapply(genome, scoreContig, scoremat=pssm, circle=T)


# For multi-contigs, usually should set circle=F
# Also requires changing subsequent code to handle lists (lapply)

# Only one contig?
np.uss.scores <- np.uss.scores[[1]]

# Check out the first few lines
head(np.uss.scores, 20)

#Check out the score dataframe
str(np.uss.scores)

#######################################
# Reduce to 1 genome length of scores #
#######################################

### Warning: Somewhat slow ###

# Get maximum score per pair of scores (w and c strands)
np.uss.scores$max    <- apply(np.uss.scores, 1, max)

# Get which strand had maximum score
np.uss.scores$strand <- max.strand(np.uss.scores)

write.csv(np.uss.scores, file="np.uss.scores.fixed.csv", quote=FALSE) #export uss scores

np.uss.scores<- read.csv("np.uss.scores.csv")


