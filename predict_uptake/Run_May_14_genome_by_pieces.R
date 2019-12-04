#######################################################################################
##  This is Program Run_May14_genome_by_pieces.R                                     ##
##  It's a model predicting uptake of DNA fragments by naturally competent bacteria. ##
##  It was created by Marcelo Mora as part of his PhD thesis work at UBC.            ##
##  Publication info:                                                                ##
##  Github location:                                                                 ##
##  This script runs the model, using settings specified in                          ##
##  Uptake_model_settings.txt and functions specified in Functions_M1_OR_2.R are use ##
##  for running model 1 or 2 and Functions_m3.R  for model 3                         ##
##  See associated ReadME file for more information.                                 ##
##  Note: this script version breaks down the genome into 100kb pieces. Then, each   ## 
##  piece runs one by one until the whole genome uptake is predicted                 ##
##  splitting the genome is useful for running the model to the entire genome        ##
##             USE THIS SCRIPT IF YOU ARE RUNNING ENTIRE GENOME                 K    ##
#######################################################################################


################################################################################
#                          load packages                                       #
################################################################################
library(dplyr)                                                                 
library(tidyverse)                                                                 
library(data.table)                                                 
library(readxl)
library(gtools)


###################################################
#                 Read in the Settings            #
###################################################

settings <- read_excel("Commons/Uptake_model_settings.xlsx")


USS_list <- settings$USS_list
folder_name <- settings$folder_name
Genome_used <- settings$Genome_used
Custom_genome_length <- settings$Custom_genome_length
Segment_start <- settings$Segment_start
Segment_end <- settings$Segment_end
Uptake_table <- settings$Uptake_table
Baseline_binding <- settings$Baseline_binding
Baseline_uptake <-  settings$Baseline_uptake
Length_of_end_overlap <- settings$Length_of_end_overlap
Fragment_distribution_file <- settings$Fragment_distribution_file
binding_prob<-settings$Binding_table
Run_ID <- paste(format(Sys.Date(), format="%B_%d_%Y"),settings$Run_ID, sep = "")
Function_file <- settings$Function_file
p_b_uss <- settings$p_b_uss


# Write the settings information to a file whose name includes the RunID:
settings_out <- data.frame(USS_list, Segment_start, Segment_end, Uptake_table, 
                           Baseline_binding, Baseline_uptake, Length_of_end_overlap, 
                           Fragment_distribution_file, Run_ID, Function_file, p_b_uss)

write.csv(settings_out, file = paste(folder_name,"settings",Run_ID, ".csv", sep = ""))

################################################################################
#                          Load the datasets                                   #
################################################################################

USS_list  <- fread(USS_list )

Uptake_table <- fread(Uptake_table)

Fragment_distribution_file <- fread(Fragment_distribution_file)

binding_prob <- fread(binding_prob)

sigmoidal_binding <- fread(binding_p) 

p_b_uss <- fread(p_b_uss)

#########################################################################################
#             Make the USS list terminally redundant                                    #
# Add the last 10 USS to the beginning of the USS list (with negative position numbers) #
# and the first 10 USS to the end of the list (with higher position numbers).           #
#                                                                                       #
#  This is done to allow uptake predictions for USSs that are close to an end.          #      
#########################################################################################


###########################
#  defines genome length  #
###########################

genome_length <- c()

if (Genome_used == "NP") {
  genome_length <- 1914386
}

if (Genome_used == "GG") {
  genome_length <- 1887046
}

if (Genome_used == "Rd") {
  genome_length <- 1831585
}

if (Genome_used == "other") {
  genome_length <- Custom_genome_length
}

# add an error if neither "NP, GG, RD or other" are chosen in Genome_used settings
if (length(genome_length) == 0) {
  print("ERROR invalid option chosen in 'Genome_used' settings. Remember, if 'other' option is used then 'Custom_genome_length' settings must be filled")
}

# Create a terminal redundant list to deal with USS close to the ends of the genome
redundant_uss_list <- rbind(tail(USS_list, Length_of_end_overlap), 
                         USS_list, head(USS_list, Length_of_end_overlap)) 

# Replace the genomic position corresponding to the position 16  of each USS 
# by with negative position numbers for the last USSs to the beginning of the USS list and with 
# higher positive numbers for the first USSs at the end of the positions at the USS list
s0 <- head(USS_list$pos, Length_of_end_overlap)                             
s1 <- tail(USS_list$pos, Length_of_end_overlap)                            
s2 <- (genome_length + s0)                                        
s3 <- -(genome_length - s1)                                       
redundant_pos <- c(s3,USS_list$pos, s2)                                
redundant_uss_list$pos <- redundant_pos                                            

# fix the number of decimals of USS score
redundant_uss_list$USS.score <- round(redundant_uss_list$USS.score, digits = 2)

USS_list <- redundant_uss_list # rename redundant USS list

#############################################################################
# load the functions to be used
#############################################################################


source(Function_file)


#  remove columns added by fread function
Fragment_distribution_file$V1 <- NULL

################################################################################ 
#     Call the model functions to predict uptake       #
################################################################################

##############################################################################
# HOW THE UPTAKE PREDICTION FUNCTIONS WORK:
# The model does a full calculation of uptake contributions all the overlapping 
# fragments only for the first focal position (in function uptake_first_pos).
# For each subsequent focal position, the model does not calculate the contributions 
# of all the overlapping fragments. Instead it take the previous position's sum and 
# modifies it by adding the contribution of the new rightmost fragment
# and subtracting the contribution of the former leftmost fragment (in function right_pe).
# It cannot efficiently save the contribution of the formerly leftmost fragment so 
# it recalculates it (in function left_pe). Function model first calls uptake_first_pos 
# to calculate the first position, and then iterates through the other positions 
# using function right_pe.
##############################################################################

genome <- c(Segment_start:Segment_end)  # pick a genomic region to be analyzed

# factor to split sequence
t<- ceiling(seq_along(genome)/100000)
# split into 20 chunks
list.genome<- split(genome, t)

for(i in 1:length(list.genome)){ 
  genome<- list.genome[[i]]
  exp <- rep(0,length(list.genome[[i]]))
  predic <- rep(0,length(list.genome[[i]]))  # create an empty vector to preallocate space
  # For each fragment size class, calculate predicted uptake for all genomic 
  # positions in vector "genome"
  for (ds in 1:length(Fragment_distribution_file$bases)) { # for each fragment size class
    # calculate predicted uptake for a genomic region using a given fragment size class
    exp.1 <- all_genome_uptake(ds) 
    exp.2 <- exp.1 * Fragment_distribution_file$frequency[ds] # multiply by frequency of that size class
    predic <- predic + exp.2
  }
  # save. This predicted uptake is not normalized yet
  write.csv(predic, file = paste(folder_name,"not_normalized_predicted_uptake",i, Run_ID, ".csv",sep = ""))
}


###################################
#   Normalize predicted uptake    #
###################################

# DO NOT MOVE THE NON-NORMALIZED UPTAKE FILES SAVED IN THE DIRECTORY BY THE PREVIOUS FUNCTION 
# OR THIS SCRIPT WON'T WORK


# creates the list of all the csv files in the directory
ldf <- list() # creates a list
listcsv <- dir(path = folder_name, pattern = paste("not_normalized_predicted_uptake",pattern = "[[:digit:]]", Run_ID, ".csv",sep = "")) 
listcsv2 <- dir(path = folder_name, pattern = paste("not_normalized_predicted_uptake",pattern = "[[:digit:]]", pattern = "[[:digit:]]", Run_ID, ".csv",sep = "")) 
listcsv3 <- c(listcsv, listcsv2)
list.name<- mixedsort(listcsv3)
# read the list of results from the folder
for (k in 1:length(list.name)){
  ldf[[k]] <- read.csv(paste(folder_name, list.name[k], sep = ""))
}

# merge all 100kb fragments results into one vector of the entire genome
data <- c()
for(i in 1:length(ldf)){
  temp<- ldf[[i]]$x
  data<- c(data,temp)
}  


# function for normalizing predicted uptake
  norm <- function(data = data) {
    h_pe <- mean(data, na.rm = FALSE)
    s_pe <- (data * 1)/h_pe
    return(s_pe)
  }
  
# normalize predicted uptake to a mean of 1
norm_exp <- norm(data = data)

# check if it was normalize correctly. Mean should be 1

mean(norm_exp)

  
# save. This predicted uptake is already normalized
write.csv(norm_exp, file = paste(folder_name,"normalized_predicted_uptake", Run_ID, ".csv",sep = ""))


