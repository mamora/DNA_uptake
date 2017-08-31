#####################################################################################
########       Calculate the distance to the closest USS                  ###########
#####################################################################################

###########################################################
######   load samples and working directory        ########         
###########################################################

# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/" # Where these files are located

# Set working directory
setwd(whereami)

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")


Uptake.ratio.np<- read.csv("./datasets/Uptake.ratio.np.csv") #read scores file 

Uptake.ratio.gg<- read.csv("./datasets/Uptake.ratio.gg.csv") #read scores file 

Uptake.uss.10.list.gg<- read.csv("./datasets/Uptake.uss.10.list.gg.csv")

Uptake.uss.10.list.np<- read.csv("./datasets/Uptake.uss.10.list.np.csv")

###########################################################################
######   calculate distance of each position to closest USS        ########         
###########################################################################


close.USS.np<- sapply(Uptake.ratio.np$pos,dist.USS, USS.genome.c = Uptake.uss.10.list.np$keypos) #calculate closest distance to USS for each position

close.USS.gg<- sapply(Uptake.ratio.gg$pos,dist.USS, USS.genome.c = Uptake.uss.10.list.gg$centralpos) #calculate closest distance to USS for each position

#add them to dataframes

Uptake.ratio.np$close.USS.np<- close.USS.np

Uptake.ratio.gg$close.USS.gg<- close.USS.gg

#save files

write.csv(Uptake.ratio.np, "./datasets/Uptake.ratio.np.csv")

write.csv(Uptake.ratio.gg, "./datasets/Uptake.ratio.gg.csv")


############################################################################