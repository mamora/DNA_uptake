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
folder.name<- "./Marcelo_paper_figures_new_order/model/Sept_24_2018/"  
#read uptake file
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")       
# list of USS10
Up.USS.np.10.list<- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) 
# fragment sizes
up15.dist1<- fread(paste(folder.name,"up15.dist_raw_den_class.csv", sep = ""))


###################################################
#                 Settings                        #
###################################################
# Settings
settings<- read.csv("./Marcelo_paper_figures_new_order/model/Sept_24_2018/Settings.csv")
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
circle.uss.list<- rbind(tail(Up.USS.np.10.list, 10), 
                        Up.USS.np.10.list, head(Up.USS.np.10.list, 9)) 

# Replace the genomic position corresponding to the pos 16 (keypos) of each uss 
# by the pretended positions if the genome was larger and a set of negative 
# positions at the beggining of the genome.                     
# This is done to facilitate calculation to find a USS in a circularized genome
s0<- head(Up.USS.np.10.list$keypos, 9) #first 9 uss                            
s1<- tail(Up.USS.np.10.list$keypos, 10) #last 9 uss                            
s2<- (length(Uptake.ratio.np$pos) + s0)                                        
s3<- -(length(Uptake.ratio.np$pos) - s1)                                       
circle.uss<- c(s3,Up.USS.np.10.list$keypos, s2)                                
circle.uss.list$keypos<- circle.uss                                            



################################################################################ 
#     Run the model funtions to calculate binding and uptake initiation        #
################################################################################
b = circle.uss.list$keypos # rename USS list as parameter "b"

# load the functions to be used
source("./Marcelo_paper_figures_new_order/model/Sept_24_2018/model_code_simple_model_v1.1.R")       

up15.dist1$bases<- as.integer(up15.dist1$bases)
up15.dist1$V1<- NULL

# make a vector with genomic positions to be used in the model
st = 1                                                                      
end = 100000
genome<- c(st:end)

# preallocate the size of the matrix to be used to store predicted uptake
f.tem<- data.frame(matrix(nrow = length(genome), ncol = length(up15.dist1$bases)))

# For each fragment size in matrix "up15.dist1", calculate predicted uptake for all genomic 
# positions in vector "genome"
for(ds in 1:length(up15.dist1$bases)){ 
tab.1<- func.model(ds)   
sum<- sapply(X = tab.1, FUN = sum.total, sim = sim)
f.tem[,ds] <- sum 
}

write.csv(f.tem, file = paste(folder.name,"simple_model_100kb_np_dist_raw_den_class.csv", sep = "")) 


