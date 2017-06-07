#####################################################################################
########       Calculate the distance to the closest USS                  ###########
#####################################################################################

###########################################################
######   load samples and working directory        ########         
###########################################################

# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/documentos phd/experiments/virtual experiments/Vex4/tutorialRjosh/Tutorial1/" # Where these files are located

# Set working directory
setwd(whereami)

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("pssmFunctions1.R")


NP.uptake.ratio<- read.csv("./datasets/NP.uptake.ratio.csv") #load uptake ratios np short fragments

PittGG.uptake.ratio<- read.csv("./datasets/PittGG.uptake.ratio.csv") #load uptake ratios PittGG short and large fragments

Up.USS.PittGG.10.list<- read.csv("./datasets/Up.USS.PittGG.10.list.csv")

Up.USS.np.10.list<- read.csv("./datasets/Up.USS.np.10.list.csv")

###########################################################################
######   calculate distance of each position to closest USS        ########         
###########################################################################


close.USS.np<- sapply(NP.uptake.ratio$pos,dist.USS, USS.genome.c = Up.USS.np.10.list$keypos) #calculate closest distance to USS for each position

close.USS.gg<- sapply(NP.uptake.ratio$pos,dist.USS, USS.genome.c = Up.USS.PittGG.10.list$keypos) #calculate closest distance to USS for each position

#add them to dataframes

NP.uptake.ratio$close.USS.np<- close.USS.np

PittGG.uptake.ratio$close.USS.gg<- close.USS.gg

#save

write.csv(NP.uptake.ratio, "./datasets/NP.uptake.ratio.csv")

write.csv(PittGG.uptake.ratio, "./datasets/PittGG.uptake.ratio.csv")


############################################################################