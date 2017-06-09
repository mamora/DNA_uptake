
###############################################################
########    Generate a list of isolated USS            ########         
###############################################################

# The objective is to generate a list of isolated USS. This list will be used to assess which position woud better represent the centre of the USS


###########################################################
######   load samples and working directory        ########         
###########################################################


# Name the path of the working directory; 
### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "~/DNA_uptake/" # Where these files are located

# Set directory
setwd(whereami)

list.files(whereami) #see files in directory

# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")

Up.USS.np.10.list<- read.csv(file="./datasets/Up.USS.np.10.list.csv")


##################################################################
#########     Generate the list of isolated uss       ############ 
##################################################################


next.USS<- sapply(Up.USS.np.10.list$keypos,dist.next.USS, USS.genome.c = Up.USS.np.10.list$keypos ) #calculate how close the USS are to next USS

Up.USS.np.10.list$next.USS<- next.USS #add to dataframe

write.csv(Up.USS.np.10.list, file="./datasets/Up.USS.np.10.list.csv", quote=FALSE) #save file

isolated.uss<- which(Up.USS.np.10.list$next.peak > 1000) # get which positions have uss with a distance higher than 1000bp from the next uss

isolated.uss.list<- Up.USS.np.10.list[isolated.uss,] # generate the list subseting the positions chosen

write.csv(isolated.uss.list, file="./datasets/isolated.uss.list.csv", quote=FALSE) #save file

###################################################################








