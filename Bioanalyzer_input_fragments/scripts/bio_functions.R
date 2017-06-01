
#############################load packages################################
library(plyr)
library(ggplot2)
library(dplyr)
library(DT)
library(cowplot)
#######################load input bioanalizer data#########################
UP13<- read.csv(file = "C:/Users/marcelo/Dropbox/uptake/Everything_else/Bioanalyzer images/input_excel/UP13.csv")
UP14<- read.csv(file = "C:/Users/marcelo/Dropbox/uptake/Everything_else/Bioanalyzer images/input_excel/UP14.csv")
UP15<- read.csv(file = "C:/Users/marcelo/Dropbox/uptake/Everything_else/Bioanalyzer images/input_excel/UP15.csv")
UP16<- read.csv(file = "C:/Users/marcelo/Dropbox/uptake/Everything_else/Bioanalyzer images/input_excel/UP16.csv")
###########################################################################


#function to add a sample column with sample name to each dataframe
add.names<- function(l.df){
for (i in 1:4){
  sample<- rep(names(l.df[i]),length(l.df[[i]]$Time))
  l.df[[i]]$sample<- sample
}
return(l.df)
  }    
  
#function to add a sample column with sample name to each dataframe
add.bases<- function(l.df){
  for (i in 1:4){
    l.df[[i]]$bases<- predict(model, newdata = l.df[[i]]$Time)  }
  return(l.df)
}    

#function to add a column with relative number of molecules
add.rel_molecules<- function(l.df){
  for (i in 1:4){
    l.df[[i]]$rel_mol<- l.df[[i]]$Value/l.df[[i]]$bases
}    
  return(l.df)
  }