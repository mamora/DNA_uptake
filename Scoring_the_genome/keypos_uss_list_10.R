#######################################################################################
#              The goal of this script os to align USSs by position 16                #
#######################################################################################

###########################################################
#####                86-028NP                      ########         
###########################################################

library(data.table)

################Load Functions################################
# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")

# Load list of uptake ratios 
Uptake.ratio.np<- fread("/datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv") #load new input samples fro UP01   
Uptake.ratio.gg<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.gg.pilon.csv") #load new input samples fro UP01   

# Load list of positions identified as USSs
Uptake.uss.10.list.np<- read.csv("./datasets/final_datasets/Up.USS.NP.10.list.pilon.csv")
Uptake.uss.10.list.gg<- read.csv("./datasets/final_datasets/Up.USS.GG.10.list.pilon.csv")


###########################################################
######   chose central position and align all USS  ########         
###########################################################

# Load the scoring model
motif   <- read.csv("~/DNA_uptake/datasets/degMotif_F", stringsAsFactors=FALSE, nrows=4, row=1) # motif model

offset<- 16  # choose the position to align USS

sitelen <- dim(motif)[2] # total number of positions in motif

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


# fix any keypos that run off the end of the linearized chromosome
tooLong <- which(Uptake.uss.10.list.np$keypos > length(Uptake.ratio.np$pos))

#########################################
# Cross-tabulate uptake ratios with USS #
#########################################

str(Uptake.ratio.np)

# uptake ratio at key position
Uptake.uss.10.list.np$keyup_small <- Uptake.ratio.np$ratio_short[match(Uptake.uss.10.list.np$keypos, Uptake.ratio.np$pos)]
Uptake.uss.10.list.np$keyup_large <- Uptake.ratio.np$ratio_long[match(Uptake.uss.10.list.np$keypos, Uptake.ratio.np$pos)]

write.csv(Uptake.uss.10.list.np, file="./datasets/final_datasets/Up.USS.NP.10.list.pilon.csv", quote=FALSE) #save file

str(Uptake.uss.10.list.np)


###########################################################
#####                PittGG                        ########         
###########################################################

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

write.csv(Uptake.uss.10.list.gg, file="./datasets/final_datasets/Up.USS.GG.10.list.pilon.csv", quote=FALSE) #save file





