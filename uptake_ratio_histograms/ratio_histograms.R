#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
#########################################################################


#############################load dataframes##############################
##WARNING MEMORY INTENSIVE###  
#YOU MIGHT HAVE TO LOAD ONE AT A TIME AND REFRESH THE SESSION OFTEN
# load dataframes
Uptake.ratio.np<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.np.csv") #load uptake ratios   
Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv")    

#########################################################################



###############################    turn wide dataframes into long dataframes #########################
######################################################################################################


Uptake.ratio.np.all<- Uptake.ratio.np[,1:8] %>% tidyr::gather(sample, ratio, 3:8) #transform from wide to long dataframe format 

str(Uptake.ratio.np.all)


str(Uptake.ratio.gg)

Uptake.ratio.gg.all<- Uptake.ratio.gg[,1:8] %>% tidyr::gather(sample, ratio, 3:8) #transform from wide to long dataframe format 



########################################    plot uptake ratios histograms      #######################################



u<- unique(Uptake.ratio.np.all$sample)



for (i in 1:length(u)){
  p<- Uptake.ratio.np.all %>% dplyr::filter(sample == u[i]) %>% ggplot() +
    geom_histogram(aes(x = ratio), binwidth = 0.05, colour = "black") +
    scale_y_continuous(limits = c(0,5e+5), expand = c(0, 0))+  
    scale_x_continuous(limits = c(0,4), breaks = seq(0 , 4, 1), expand = c(0, 0))+ 
    labs(x = "uptake ratios") +
    ggtitle(paste("Histogram of uptake ratios from ", u[i],"uptake donor DNA", sep=" ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/uptake ratios/uptake_ratios_histograms/histogram","ratios","4", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 500, units = "px")
  print(p)
  dev.off()
}  


u<- unique(Uptake.ratio.gg.all$sample)


for (i in 1:length(u)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u[i]) %>% ggplot() +
    geom_histogram(aes(x = ratio), binwidth = 0.05, colour = "black") +
    scale_y_continuous(limits = c(0,5e+5), expand = c(0, 0))+  
    scale_x_continuous(limits = c(0,4), breaks = seq(0 , 4, 1), expand = c(0, 0))+ 
    labs(x = "uptake ratios") +
    ggtitle(paste("Histogram of uptake ratios from ", u[i],"uptake donor DNA", sep=" ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/uptake ratios/uptake_ratios_histograms/histogram","ratios", u[i],"4", "tiff", sep=".")
  tiff(file_name, width = 1000, height = 500, units = "px")
  print(p)
  dev.off()
}




for (i in 1:length(u)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u[i]) %>% ggplot() +
    geom_histogram(aes(x = ratio), binwidth = 0.05, colour = "black") +
    scale_y_continuous(limits = c(0,1e+4), expand = c(0, 0))+  
    scale_x_continuous(limits = c(4,20), breaks = seq(4 , 20, 2), expand = c(0, 0))+ 
    labs(x = "uptake ratios") +
    ggtitle(paste("Histogram of uptake ratios from ", u[i],"uptake donor DNA", sep=" ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/uptake ratios/uptake_ratios_histograms/histogram","ratios", u[i],"3", "tiff", sep=".")
  tiff(file_name, width = 1000, height = 500, units = "px")
  print(p)
  dev.off()
}