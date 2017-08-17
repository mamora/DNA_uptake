
#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)

# load dataframes (OPTIONAL CODE)
Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv") #load new input samples fro UP01   
Uptake.ratio.np<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.np.csv") #load new input samples fro UP01   




#####################################################################
##############################Functions##############################
#####################################################################


synthetic<- function(block, PittGG.block, np.block){

  #find positions with "n" and without "n"
gaps<- which(np.block[[1]] == "-")
nongaps<- which(np.block[[1]] != "-") 

#get size of the colinear block
np.block.size<- length(np.block[[1]])

pos <- np.block #renamed

p<- c(block$left[2]:block$right[2]) # make a vector with genomic positions of the block 

positions.np<- c(1:length(pos[[1]])) #make a vector of the same size as the block

positions.np[nongaps]<- p # fill the position vector with genomic positions for positions without gaps

g<- rep(0,length(gaps))  # create a vector of gaps of value 0

positions.np[gaps]<- g  # fill the positions with value of 0 for gaps

############################################################
####################Do the same for PittGG##################
############################################################

#find positions with "n"
gaps<- which(PittGG.block [[1]] == "-")
nongaps<- which(PittGG.block [[1]] != "-") 
GG.block.size<- length(PittGG.block [[1]])
pos <- PittGG.block 
p<- c(block$left[1]:block$right[1]) # make a vector with genomic positions of the block 
positions.GG<- c(1:length(pos[[1]])) # make a vector of the same size as the block
positions.GG[nongaps]<- p #fill the position vector with genomic positions for positions without gaps
g<- rep(0,length(gaps))  # create a vector of gaps of value 0
positions.GG[gaps]<- g  ## fill the positions with value of 0 for gaps


###########################################################
#################create matrix#############################
###########################################################

m<- matrix(data = NA, nrow = 4, ncol = length(PittGG.block[[1]]))

m[1,]<- np.block[[1]]
m[2,]<- PittGG.block[[1]]
m[3,]<- positions.np
m[4,]<- positions.GG

return(m)

}


syn_to_data<- function(m){

#get positions without gaps in NP
pos.np<- m["NP_positions",(which(m["NP_sequence",] != 6))]
#get positions without gaps in PittGG
pos.gg<- m["PittGG_positions",(which(m["PittGG_sequence",] != 6))]

#get size of colinear block
size<- dim(m)[2] 

#create a vector of NA for each position --> later will be this will be gaps
np.up<- rep(NA, size)

# get uptake ratios of non-gaps positions for NP
np.ratio.nongaps.s<- Uptake.ratio.np$ratio_short[pos.np]

# get uptake ratios of non-gaps positions for PittGG
GG.ratio.nongaps.s<- Uptake.ratio.gg$ratio_short[pos.gg]

# get uptake ratios of non-gaps positions for NP
np.ratio.nongaps.l<- Uptake.ratio.np$ratio_long[pos.np]

# get uptake ratios of non-gaps positions for PittGG
GG.ratio.nongaps.l<- Uptake.ratio.gg$ratio_long[pos.gg]



num<- c(1:size)

#create a dataframe to make an uptake figure so far with only NA instead of ratios
block.synthenic<- data.frame(relative_pos = num, NP_uptake.short = np.up , PittGG_uptake.short = np.up, NP_uptake.long = np.up , PittGG_uptake.long = np.up)


#replace NA's with uptake ratios for non-gaps positions

block.synthenic$NP_uptake.short[(which(m["NP_sequence",] != 6))]<- np.ratio.nongaps.s

block.synthenic$PittGG_uptake.short[(which(m["PittGG_sequence",] != 6))]<- GG.ratio.nongaps.s

block.synthenic$NP_uptake.long[(which(m["NP_sequence",] != 6))]<- np.ratio.nongaps.l

block.synthenic$PittGG_uptake.long[(which(m["PittGG_sequence",] != 6))]<- GG.ratio.nongaps.l



return(block.synthenic)

}

