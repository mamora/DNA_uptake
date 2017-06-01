#####################################################################################
########       align USSs positions by highest uptake position            ###########
#####################################################################################

#This analysis is the continuation from the centre_USS script


###########################################################
#####                86-028NP                      ########         
###########################################################

###########################################################
######   load samples and working directory        ########         
###########################################################

# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/documentos phd/experiments/virtual experiments/Vex4/tutorialRjosh/Tutorial1/" # Where these files are located

# Set working directory
setwd(whereami)

list.files(whereami) #see files in directory

# All this dataframes are available by request. Contact redfield@zoology.ubc.ca

################Load Functions################################

source("C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/peak finder/pssmFunctions1.R")

##########Load list of uptake ratios per genomic position#########################

uptake<- read.csv("C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Scoring_USS/Output/uptake_recal.csv") #save file


###########Load list of positions identified as USS##############################


Uptake.uss.10.list<- read.csv("C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Scoring_USS/Output/Uptake.uss.list.cut.10.csv")


###########################################################
######   chose central position and align all USS  ########         
###########################################################


offset<- 17  # set offset positions in the motif which has the highest uptake of the uss peak

#######################
# Index USS positions #
#######################

# index site orientations
f.ind <- which(Uptake.uss.10.list$strand == "w")
r.ind <- which(Uptake.uss.10.list$strand == "c")


# offset to the key position
Uptake.uss.10.list$keypos <- 0 # initialize
Uptake.uss.10.list$keypos[f.ind] <- Uptake.uss.10.list$USS.pos[f.ind] + offset - 1 # forwards
Uptake.uss.10.list$keypos[r.ind] <- Uptake.uss.10.list$USS.pos[r.ind] + sitelen - offset # reverses

str(Uptake.uss.10.list)


# fix any keypos that run off the end of the linearized chromosome
tooLong <- which(Uptake.uss.10.list$keypos > length(uptake$pos))

#########################################
# Cross-tabulate uptake ratios with USS #
#########################################

# flag at key position
Uptake.uss.10.list$keyflag_small <- uptake$flag_small[match(Uptake.uss.10.list$keypos, uptake$pos)]

Uptake.uss.10.list$keyflag_large <- uptake$flag_large[match(Uptake.uss.10.list$keypos, uptake$pos)]

# uptake ratio at key position
Uptake.uss.10.list$keyup_small <- uptake$ratio_small[match(Uptake.uss.10.list$keypos, uptake$pos)]

Uptake.uss.10.list$keyup_large <- uptake$ratio_large[match(Uptake.uss.10.list$keypos, uptake$pos)]


# need a check here to ensure the liftover is correct, but it should be...


write.csv(Uptake.uss.10.list, file="C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Scoring_USS/Output/Uptake.uss.list.cut.10.csv", quote=FALSE) #save file




###########################################################
#####                PittGG                        ########         
###########################################################

###########################################################
######   load samples and working directory        ########         
###########################################################

Up.PittGG.10.USS.list <- read.csv("C:/Users/marcelo/Dropbox/uptake/Everything_else/Marcelo scripts and files/PittGG_maps/Up.PittGG.10.USS.list.csv")


uss<- Up.PittGG.10.USS.list #just rename the dataframe


# Load the scoring model
motif   <- read.csv("C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Commons/uptake.pssm.csv", stringsAsFactors=FALSE, nrows=4, row=1) # motif model

sitelen <- dim(motif)[2] # positions in motif

offset  <- 17  # position lovated in the middle of the uptake peak


#######################
# Index USS positions #
#######################

# index site orientations
f.ind <- which(uss$strand == "w")
r.ind <- which(uss$strand == "c")

# offset to the key position
uss$centralpos <- 0 # initialize
uss$centralpos[f.ind] <- uss$USS.pos[f.ind] + offset - 1 # forwards
uss$centralpos[r.ind] <- uss$USS.pos[r.ind] + sitelen - offset # reverses

str(uss)


# Reorder based on centralpos 
uss <- uss[order(uss$centralpos), ]


# uptake ratio at key position for small and large fragments 
uss$keyup_small <- PittGG.uptake.ratio$ratio_small[match(uss$centralpos, PittGG.uptake.ratio$pos)]

uss$keyup_large <- PittGG.uptake.ratio$ratio_long[match(uss$centralpos, PittGG.uptake.ratio$pos)]


write.csv(uss, file="C:/Users/marcelo/Dropbox/uptake/Everything_else/Marcelo scripts and files/PittGG_maps/Up.PittGG.10.USS.list.csv", quote=FALSE) #save file


