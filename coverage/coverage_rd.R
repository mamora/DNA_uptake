
#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
#########################################################################



#############################load dataframes##############################
rd<-  fread("~/DNA_uptake/datasets/RR722.rd.pilon.all.depth.bed") 
###############################################################################

###############################trim an save rd depth file ############################## 
str(rd)
unique(rd$V4)
rd.all<- rd %>% filter(V4 == "gi|16271976|ref|NC_000907.1|pilon") #extract only rd 
rd.all<- rd.all[,7:8] #extract only useful columns
colnames(rd.all)<- c("pos","depth") # add headers
str(rd.all)
########################################################################################



##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
rd.all.mapped<- 10168072 # total mapped reads from summary table
rd.n<- (rd.all$depth *1e+6)/up16.all.mapped
#########################################################################

