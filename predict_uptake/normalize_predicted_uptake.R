######################################################################
##      Normalize predicted DNA uptake in NP small frag data     #####
######################################################################

###########################################################
## Warning: Please open Uptake_summer2017.Rproj first #####
###########################################################

##############################################################
######   1. load samples and working directory        ########         
##############################################################

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)

                                                               
folder.name <- "./Marcelo_paper_figures_new_order/model/Oct_16_2018_2/"  
#read uptake file
Uptake.ratio.np <- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")      
# load USS scores
USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")
# list of USS10
Up.USS.np.10.list <- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) 

###################################################################################################
#              creates the list of all the csv files in the directory.                            #
#        Each files is the non-normalized predicted uptake for a 100kb chunck                     #
###################################################################################################
ldf <- list() # creates a list                                                                    #
listcsv <- dir(path = folder.name,                                                                #
               pattern = "model_sigmodial_100kbnp_class_sizessmall[[:digit:]]*.csv")              #
list.name<- mixedsort(listcsv)                                                                    #
# read the list of results from the folder                                                        #
for (k in 1:length(list.name)){                                                                   #
  ldf[[k]] <- read.csv(paste(folder.name, list.name[k], sep = ""))                                #
}                                                                                                 #
# merge all 100kb fragments results into one vector of the entire genome                          #
up15.res<- c()                                                                                    #
for(i in 1:length(ldf)){                                                                          #
temp<- ldf[[i]]$x                                                                                 #
up15.res<- c(up15.res,temp)                                                                       #
}                                                                                                 #
###################################################################################################

# read fragment size classes
up15.dist1.class.sizes <- fread(paste(folder.name,"up15.dist_size.class.csv", sep = ""))

###################################
#   Normalize predicted uptake    #
###################################
# function for normalizing predicted uptake
norm<- function (data = data){
  h_pe<- mean(data, na.rm = FALSE)
  s_pe<- (data * 1)/h_pe
  return(s_pe)
}


# normalize predicted uptake to a mean of 1
s1.class<- norm(data = up15.res)
# normalize predicted uptake to a mean of 1
ob<- Uptake.ratio.np$ratio_short

# check if it was normalize correctly. Mean should be 1
mean(s1.class)
mean(ob)


# put observed and predicted uptake in a dataframe
real.data<- data.table(pos = Uptake.ratio.np$pos, observed_uptake = ob,
                       expected_uptake = s1.class)

# save model results
write.csv(real.data, 
          file = paste(folder.name, "small_predicted_sigmodial_model_corrected_10.csv", sep = ""))

