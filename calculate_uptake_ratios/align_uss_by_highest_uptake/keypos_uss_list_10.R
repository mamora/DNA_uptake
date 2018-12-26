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
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/" # Where these files are located

# Set working directory
setwd(whereami)

list.files(whereami) #see files in directory

# All this dataframes are available by request. Contact redfield@zoology.ubc.ca

################Load Functions################################
# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("~/DNA_uptake/helper_functions/pssmFunctions1.R")

##########Load list of uptake ratios per genomic position#########################
# load dataframes
Uptake.ratio.np<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.np.csv") #load new input samples fro UP01   
# load dataframes
Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv") #load new input samples fro UP01   

###########Load list of positions identified as USS##############################


Uptake.uss.10.list.np<- read.csv("./datasets/Up.USS.np.10.list.csv")

Uptake.uss.10.list.gg<- read.csv("./datasets/final_datasets/USS_scores/uptake_model/Up.USS.PittGG.9.5.list.pilon.csv")


###########################################################
######   chose central position and align all USS  ########         
###########################################################

# Load the scoring model
motif   <- read.csv("~/DNA_uptake/datasets/uptake.pssm.csv", stringsAsFactors=FALSE, nrows=4, row=1) # motif model


offset<- 16  # set offset positions in the motif which has the highest uptake of the uss peak

sitelen <- dim(motif)[2] # positions in motif

#######################
# Index USS positions #
#######################

# index site orientations
f.ind <- which(Uptake.uss.10.list.np$strand == "w")
r.ind <- which(Uptake.uss.10.list.np$strand == "c")


# offset to the key position
Uptake.uss.10.list.np$keypos <- 0 # initialize
Uptake.uss.10.list.np$keypos[f.ind] <- Uptake.uss.10.list.np$USS.pos[f.ind] + offset - 1 # forwards
Uptake.uss.10.list.np$keypos[r.ind] <- Uptake.uss.10.list.np$USS.pos[r.ind] + sitelen - offset # reverses

str(Uptake.uss.10.list.np)


# fix any keypos that run off the end of the linearized chromosome
tooLong <- which(Uptake.uss.10.list.np$keypos > length(Uptake.ratio.np$pos))

#########################################
# Cross-tabulate uptake ratios with USS #
#########################################

str(Uptake.ratio.np)

# flag at key position
Uptake.uss.10.list.np$flag_small <- Uptake.ratio.np$flag_small[match(Uptake.uss.10.list.np$keypos, Uptake.ratio.np$pos)]

Uptake.uss.10.list.np$flag_large <- Uptake.ratio.np$flag_large[match(Uptake.uss.10.list.np$keypos, Uptake.ratio.np$pos)]

# uptake ratio at key position
Uptake.uss.10.list.np$keyup_small <- Uptake.ratio.np$ratio_short[match(Uptake.uss.10.list.np$keypos, Uptake.ratio.np$pos)]

Uptake.uss.10.list.np$keyup_large <- Uptake.ratio.np$ratio_long[match(Uptake.uss.10.list.np$keypos, Uptake.ratio.np$pos)]


# need a check here to ensure the liftover is correct, but it should be...


write.csv(Uptake.uss.10.list.np, file="./datasets/Uptake.uss.10.list.np.csv", quote=FALSE) #save file

str(Uptake.uss.10.list.np)


###########################################################
#####                PittGG                        ########         
###########################################################

###########################################################
######   load samples and working directory        ########         
###########################################################

str(Uptake.uss.10.list.gg)


#######################
# Index USS positions #
#######################

# index site orientations
f.ind <- which(Uptake.uss.10.list.gg$strand == "w")
r.ind <- which(Uptake.uss.10.list.gg$strand == "c")

# offset to the key position
Uptake.uss.10.list.gg$centralpos <- 0 # initialize
Uptake.uss.10.list.gg$centralpos[f.ind] <- Uptake.uss.10.list.gg$USS.pos[f.ind] + offset - 1 # forwards
Uptake.uss.10.list.gg$centralpos[r.ind] <- Uptake.uss.10.list.gg$USS.pos[r.ind] + sitelen - offset # reverses

str(Uptake.uss.10.list.gg)


# Reorder based on centralpos 
Uptake.uss.10.list.gg <- Uptake.uss.10.list.gg[order(Uptake.uss.10.list.gg$centralpos), ]


# uptake ratio at key position for small and large fragments 
Uptake.uss.10.list.gg$keyup_small <- Uptake.ratio.gg$ratio_short[match(Uptake.uss.10.list.gg$centralpos, Uptake.ratio.gg$pos)]

Uptake.uss.10.list.gg$keyup_large <- Uptake.ratio.gg$ratio_long[match(Uptake.uss.10.list.gg$centralpos, Uptake.ratio.gg$pos)]


write.csv(Uptake.uss.10.list.gg, file="datasets/Uptake.uss.10.list.gg.csv", quote=FALSE) #save file





