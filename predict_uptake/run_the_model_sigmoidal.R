######################################################################
## create a model that predicts DNA uptake in NPlarge frag data #####
######################################################################

# the objective of this script is to run the model to predict uptake


################################################################################
#                          load packages                                       #
################################################################################
library(dplyr)                                                                 
library(tidyr)                                                                 
library(data.table)                                                            
################################################################################
#                          load datasets                                       #
################################################################################
# uptake ratios                                                                
folder.name <- "./Marcelo_paper_figures_new_order/model/Oct_16_2018_2/"  
# read uptake file
Uptake.ratio.np <- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")      
# list of USS10
Up.USS.np.10.list <- fread("./datasets/final_datasets/Uptake.uss.10.list.np.csv") 
# fragment sizes
#up15.dist1 <- fread(paste(folder.name,"up15.dist_all_frag_sizes.csv", sep = ""))
up15.dist1 <- fread(paste(folder.name,"up15.dist_size.class.csv", sep = ""))

###################################################
#                 Settings                        #
###################################################
# Settings
settings <- read.csv("./Marcelo_paper_figures_new_order/model/Oct_16_2018_2/Settings.csv")
# set the background binding and uptake probabilities          
b0 = settings$b0
u0 = settings$u0                                  

################################################################################
#                    circularize the USS list                                  #
# Add the last 10 USS to the beginning of the uss list and the first 9 USS to  #
# the end of the list.                                                         #
#                                                                              #
# Then replace the genomic position corresponding to the pos 16 (keypos) of    #
# each uss list by pretended positions if the genome was larger. This is done  #
# to facilitate futher calculations                                            #
################################################################################

# Add the 10 last uss and the first 9 uss to the uss list                     
circle.uss.list <- rbind(tail(Up.USS.np.10.list, 10), 
                        Up.USS.np.10.list, head(Up.USS.np.10.list, 9)) 

# Replace the genomic position corresponding to the pos 16 (keypos) of each uss 
# by the pretended positions if the genome was larger and a set of negative 
# positions at the begining of the genome.                     
# This is done to facilitate calculation to find a USS in a circularized genome
s0 <- head(Up.USS.np.10.list$keypos, 9) #first 9 uss                            
s1 <- tail(Up.USS.np.10.list$keypos, 10) #last 9 uss                            
s2 <- (length(Uptake.ratio.np$pos) + s0)                                        
s3 <- -(length(Uptake.ratio.np$pos) - s1)                                       
circle.uss <- c(s3,Up.USS.np.10.list$keypos, s2)                                
circle.uss.list$keypos <- circle.uss                                            



################################################################################ 
#     Run the model funtions to calculate binding and uptake initiation        #
################################################################################
b = circle.uss.list$keypos # rename USS list as parameter "b"

circle.uss.list$USS.score <- round(circle.uss.list$USS.score, digits = 1)

# load the functions to be used
source("./Marcelo_paper_figures_new_order/model/Oct_16_2018_2/model_code_sigmodial_model_v2.R")

up15.dist1$V1 <- NULL
# rename columns
names(up15.dist1)[1] <-"bases"
names(up15.dist1)[2] <-"ave.density"

#############################################################################
#   split the genome in chucks to save results as 100kb individual files    #
#############################################################################
# make a vector with genomic positions to be used in the model              #
st = 1                                                                      #
end = 1914386                                                               #
genome <- c(st:end)                                                         #
                                                                            #
# split the genome into a 100kb chuck elements of a list                    # 
t<- ceiling(seq_along(genome)/100000)                                       #
list.genome<- split(genome, t)                                              #
#############################################################################


for(i in 1:length(list.genome)){ 
  genome<- list.genome[[i]]  # pick a 100kb region chunck
  exp <- rep(0,length(list.genome[[i]]))  # create an empty vector to preallocate space
  # For each fragment size in matrix "up15.dist1", calculate predicted uptake for all genomic 
  # positions in vector "genome"
  for (ds in 1:length(up15.dist1$bases)) { # for each fragment size
    # calculate predicted uptake for a given 100kb chunck of the genome
    exp.1 <- func.model(ds) 
    exp.2 <- exp.1 * up15.dist1$ave.density[ds] # multiply 
    exp <- exp + exp.2
  }
  # save as 100kb chucks individual files. This predicted uptake is not normalized yet
write.csv(exp, file = paste(folder.name,"model_sigmodial_100kb",
                            "np","_class_sizes","small",i,".csv",sep = ""))
}

