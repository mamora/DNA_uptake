#####################################################################################
########       Scoring Rd genome with Uptake-bias motif model       ###########
#####################################################################################


###########################################################
######   load samples and working directory        ########         
###########################################################


# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/" # Where these files are located

# Name the sites file, and genome reference file
fasta    <- "datasets/rd.pilon.fasta" # Reference fasta to the NP genome

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
rd.uptake.scores <- lapply(genome, scoreContig, scoremat=uptake.pssm, circle=T)


# For multi-contigs, usually should set circle=F
# Also requires changing subsequent code to handle lists (lapply)

# Only one contig?
rd.uptake.scores <- rd.uptake.scores[[1]]

# Check out the first few lines
head(rd.uptake.scores, 50)

#Check out the score dataframe
str(rd.uptake.scores)

#######################################
# Reduce to 1 genome length of scores #
#######################################

### Warning: Somewhat slow ###

# Get maximum score per pair of scores (w and c strands)
rd.uptake.scores$max  <- apply(rd.uptake.scores, 1, max)

# Get which strand had maximum score
rd.uptake.scores$strand <- max.strand(rd.uptake.scores)

View(rd.uptake.scores) #check the score dataframe

write.csv(rd.uptake.scores, file="datasets/rd.uptake.scores.csv", quote=FALSE) #export uss scores to your working directory

rd.uptake.scores<- read.csv("rd.uptake.scores.csv") #read the file in case I am starting all over


rd.uptake.scores[which.min(rd.uptake.scores$w),] #check lowest score

rd.uptake.scores[which.max(rd.uptake.scores$w),] #check lowest score


rd.uptake.scores[which.min(rd.uptake.scores$c),] #check highest score

rd.uptake.scores[which.max(rd.uptake.scores$c),] #check highest score

#########################################################################################
#####################      Plot histograms of uptake scores         #####################
#########################################################################################

p <- ggplot() +
  geom_histogram(aes(x = w), binwidth = 0.2, colour = "black", data = rd.uptake.scores) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.2e+5), expand = c(0, 0))+
  labs(x = "scores Watson strand") +
  ggtitle("Histogram of uptake scores for Rd genome Watson strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","RD", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = w), binwidth = 0.2, colour = "black", data = rd.uptake.scores) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 2e+3), expand = c(0, 0))+
  labs(x = "scores Watson strand") +
  ggtitle("Histogram of uptake scores for Rd genome Watson strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","Rd", "zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = c), binwidth = 0.2, colour = "black", data = rd.uptake.scores) +
  scale_x_continuous(limits = c(-2, 13), breaks = seq(-2 , 13, 1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.2e+5), expand = c(0, 0))+
  labs(x = "scores Crick strand") +
  ggtitle("Histogram of uptake scores for Rd genome Crick strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","Rd","c", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = c), binwidth = 0.2, colour = "black", data = rd.uptake.scores) +
  scale_x_continuous(limits = c(8, 13), breaks = seq(8 , 13, 0.5))+
  scale_y_continuous(limits = c(0, 2e+3), expand = c(0, 0))+
  labs(x = "scores Crick strand") +
  ggtitle("Histogram of uptake scores for Rd genome Crick strand") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/scores/scores","Rd","c", "zoom", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()

