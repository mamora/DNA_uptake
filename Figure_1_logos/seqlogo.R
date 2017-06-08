
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
library(RColorBrewer)

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
pfm<- pfm[,4:34]

#make the PPM or PFM matrix used for the logo
p<- makePWM(pfm)



###########################################################
######   Remake the Uptake motif logo              ########         
###########################################################


setwd("C:/Users/marcelo/Dropbox/ Uptake ms/DegenerateStuff/")

#### Contains a bunch of functions (probably hacky)
# Main one is SeqLogo2, which is my hack to get the right object 
load("Outputs/functionsMarch11.R")

inpPCM <- read.csv("Outputs/inppcm.csv", row.names=1)
perPCM <- read.csv("Outputs/perpcm.csv", row.names=1)

# This gets the consensus of the right 31 bases and makes a numeric string too.
consensus <- "ATGCCAAAGTGCGGTTAATTTTTACAGTATTTTTGGGTTCGA"
consensus <- substring(consensus, 6, 36)
connum    <- as.numeric(factor(unlist(strsplit(consensus, ""))))
degen     <- 0.24 # This is the predicted dengeneracy of the input per position
sampn     <- 1e7  # This is the subsample size to make both datasets the same size.

# inpcons <- "AAAGTGCGGTTAATTTTTACAGTATTTTTGG"
# connum
# as.numeric(factor(unlist(strsplit(inpcons, ""))))

# subset looks awesome
# so does extract

# Expected Input... hacky as can be.  2010.
expPCM           <- matrix(degen/3, length(connum), 4)
rownames(expPCM) <- rownames(inpPCM)
colnames(expPCM) <- colnames(inpPCM)
for(i in 1:31){
  expPCM[i, connum[i]] <- 0.76
}
expPCM <- sampn * expPCM

# WebLogos

# This is how the input differs from expected
z <- inpPCM/expPCM / rowSums(inpPCM/expPCM)
seqLogo2(t(z))
sum(pwm2ic2(t(z))) # bits inp vs exp


# This is the final derived motif,
# I.e. this is the "uptake motif" based on the degenerate uptake experiment.
# Essentially, this says what the motif would look like IF the bases had been equiprobable at each position,
# It accounts for the position-specific degeneracy.
y <- perPCM/inpPCM / rowSums(perPCM/inpPCM)
seqLogo2(t(y))
sum(pwm2ic2(t(y))) # bits per vs inp

# Thisis how the motif would look, if we'd just assumed the input degeneracy was exact
x <- perPCM/expPCM / rowSums(perPCM/expPCM)


###########################################################
######              Plot both logos                ########         
###########################################################

png(filename="C:/Users/marcelo/Documents/DNA_uptake/Figure_1_logos/genomic_motifs.png", width = 900, height = 500, units = "px")
seqLogo(p)
dev.off()


png(filename="C:/Users/marcelo/Documents/DNA_uptake/Figure_1_logos/uptake_motif.png", width = 900, height = 500, units = "px")
seqLogo2(t(x))
dev.off()











  