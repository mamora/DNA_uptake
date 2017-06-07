
###########################################################
######   Remake the Genomic motif logo             ########         
###########################################################

# The objective is to make the genomic motif logo based on Maughan et al. [2010]


###########################################################
######   load samples and working directory        ########         
###########################################################


# Name the path of the working directory; 
### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/datasets/" # Where these files are located


# Set directory
setwd(whereami)

# Name the sites file, and genome reference file
sites    <- "Hin2206sites.txt" # From Rosie's Gibbs sampling of rd.fa

list.files(whereami) #see files in directory

library(seqLogo)

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("../functions/pssmFunctions1.R")


#################################################################
##########Remake the Genomic USS motif Logo###################### 
#################################################################

#############
# Make PSSM #
#############

sitemat <- sites2matrix(sites)    # set of sites
pcm     <- matrix2pcm(sitemat)    # count matrix
pfm     <- pcm2pfm(pcm, pseudo=0) # frequency matrix

# take only positions I used for the analysis
pfm<- pfm[,3:33]

#make the PPM or PFM matrix used for the logo
p<- makePWM(pfm)

#make the logo
seqLogo(p)






  