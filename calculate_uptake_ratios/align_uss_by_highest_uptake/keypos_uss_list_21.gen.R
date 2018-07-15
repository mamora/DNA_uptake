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
Uptake.ratio.np<- fread("./datasets/new_norm_datasets/Uptake.ratio.np.corrected.csv") #read uptake file 
# load dataframes
Uptake.ratio.gg<- fread("./datasets/new_norm_datasets/Uptake.ratio.gg.corrected.csv") #load new input samples fro UP01   

###########Load list of positions identified as USS##############################


Gen.USS.NP.21.list<- read.csv("datasets/ Gen.USS.NP.21.list.csv")

Gen.USS.PittGG.21.list<- read.csv("datasets/ Gen.USS.PittGG.21.list.csv")


###########################################################
######   chose central position and align all USS  ########         
###########################################################

# Load the scoring model
motif   <- read.csv("datasets/pssm.genomic.csv", stringsAsFactors=FALSE, nrows=4, row=1) # motif model

offset<- 16  # set offset positions in the motif which has the highest uptake of the uss peak

sitelen <- dim(motif)[2] # positions in motif

#######################
# Index USS positions #
#######################

# index site orientations
f.ind <- which(Gen.USS.NP.21.list$strand == "w")
r.ind <- which(Gen.USS.NP.21.list$strand == "c")


# offset to the key position
Gen.USS.NP.21.list$keypos <- 0 # initialize
Gen.USS.NP.21.list$keypos[f.ind] <- Gen.USS.NP.21.list$USS.pos[f.ind] + offset - 1 # forwards
Gen.USS.NP.21.list$keypos[r.ind] <- Gen.USS.NP.21.list$USS.pos[r.ind] + sitelen - offset # reverses

str(Gen.USS.NP.21.list)


# fix any keypos that run off the end of the linearized chromosome
tooLong <- which(Gen.USS.NP.21.list$keypos > length(Uptake.ratio.np$pos))

#########################################
# Cross-tabulate uptake ratios with USS #
#########################################

# uptake ratio at key position
Gen.USS.NP.21.list$keyup_small <- Uptake.ratio.np$ratio_short[match(Gen.USS.NP.21.list$keypos, Uptake.ratio.np$pos)]

Gen.USS.NP.21.list$keyup_large <- Uptake.ratio.np$ratio_long[match(Gen.USS.NP.21.list$keypos, Uptake.ratio.np$pos)]


# need a check here to ensure the liftover is correct, but it should be...


write.csv(Gen.USS.NP.21.list, file="./datasets/Gen.USS.NP.21.list.csv", quote=FALSE) #save file

str(Gen.USS.NP.21.list)


###########################################################
#####                PittGG                        ########         
###########################################################
###########Load list of positions identified as USS##############################

Gen.USS.PittGG.21.list<- read.csv("datasets/ Gen.USS.PittGG.21.list.csv")


###########################################################
######   chose central position and align all USS  ########         
###########################################################

# Load the scoring model
motif   <- read.csv("datasets/pssm.genomic.csv", stringsAsFactors=FALSE, nrows=4, row=1) # motif model

offset<- 16  # set offset positions in the motif which has the highest uptake of the uss peak

sitelen <- dim(motif)[2] # positions in motif

#######################
# Index USS positions #
#######################

# index site orientations
f.ind <- which(Gen.USS.PittGG.21.list$strand == "w")
r.ind <- which(Gen.USS.PittGG.21.list$strand == "c")


# offset to the key position
Gen.USS.PittGG.21.list$keypos <- 0 # initialize
Gen.USS.PittGG.21.list$keypos[f.ind] <- Gen.USS.PittGG.21.list$USS.pos[f.ind] + offset - 1 # forwards
Gen.USS.PittGG.21.list$keypos[r.ind] <- Gen.USS.PittGG.21.list$USS.pos[r.ind] + sitelen - offset # reverses

str(Gen.USS.PittGG.21.list)


# fix any keypos that run off the end of the linearized chromosome
tooLong <- which(Gen.USS.PittGG.21.list$keypos > length(Uptake.ratio.gg$pos))

#########################################
# Cross-tabulate uptake ratios with USS #
#########################################

# uptake ratio at key position
Gen.USS.PittGG.21.list$keyup_small <- Uptake.ratio.gg$ratio_short[match(Gen.USS.PittGG.21.list$keypos, Uptake.ratio.gg$pos)]

Gen.USS.PittGG.21.list$keyup_large <- Uptake.ratio.gg$ratio_long[match(Gen.USS.PittGG.21.list$keypos, Uptake.ratio.gg$pos)]


# need a check here to ensure the liftover is correct, but it should be...


write.csv(Gen.USS.PittGG.21.list, file="./datasets/Gen.USS.PittGG.21.list.csv", quote=FALSE) #save file

str(Gen.USS.PittGG.21.list)

#######################################################################
#########   compare USS list from uptake and genomic models  ##########
#######################################################################

Up.USS.np.9.5.list<- fread("./datasets/Up.USS.np.9.5.list.csv") #save file

Gen.USS.NP.21.list$keypos



c<- 

c2<- 

(which(c %in% c2))


length(which(Gen.USS.NP.21.list$keypos %in% Up.USS.np.9.5.list$keypos))

length(which(Up.USS.np.9.5.list$keypos %in% Gen.USS.NP.21.list$keypos))

match.up<- (which(Up.USS.np.9.5.list$keypos %in% Gen.USS.NP.21.list$keypos))


non.match.upt<- Up.USS.np.9.5.list[-match.up]

summary(non.match.upt$keyup_small)

summary(Up.USS.np.9.5.list$keyup_small)


1586/1607

1586/2248

