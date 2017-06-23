
# load uptake ratio dataframes for PittGG and NP

PittGG.uptake.ratio<- read.csv("C:/Users/marcelo/Dropbox/uptake/Everything_else/Marcelo scripts and files/PittGG_maps/PittGG.uptake.ratio.csv") #save file

np.uptake<- read.csv("C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Scoring_USS/Output/uptake_recal.csv") #save file


uptake.kb.6<- read.csv("C:/Users/marcelo/Dropbox/uptake/Final_resources/R_scripts/Scoring_USS/Output/uptake.kb.6.csv") #save file


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

p<- c(block$left[1]:block$right[1]) # make a vector with genomic positions of the block 

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
p<- c(block$left[2]:block$right[2]) # make a vector with genomic positions of the block 
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
np.ratio.nongaps<- np.uptake$ratio[pos.np]

# get uptake ratios of non-gaps positions for PittGG
GG.ratio.nongaps<- PittGG.uptake.ratio$ratio_small[pos.gg]

num<- c(1:size)

#create a dataframe fot for making an uptake figure so far with only NA instead of ratios
block.synthenic<- data.frame(relative_pos = num, NP_uptake = np.up , PittGG_uptake = np.up)


#replace NA's with uptake ratios for non-gaps positions
block.synthenic$NP_uptake[(which(m["NP_sequence",] != 6))]<- np.ratio.nongaps

block.synthenic$PittGG_uptake[(which(m["PittGG_sequence",] != 6))]<- GG.ratio.nongaps

return(block.synthenic)

}

