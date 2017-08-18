
#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
#########################################################################


#############################load dataframes##############################
##WARNING MEMORY INTENSIVE###  
#YOU MIGHT HAVE TO LOAD ONE AT A TIME AND REFRESH THE SESSION OFTEN
# load dataframes
Uptake.ratio.np<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.np.csv") #load new input samples fro UP01   
Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv") #load new input samples fro UP01   

#########################################################################



###############################    turn wide dataframes into long dataframes #########################
######################################################################################################


Uptake.ratio.np.all<- Uptake.ratio.np[,1:8] %>% tidyr::gather(sample, ratio, 3:8) #transform from wide to long dataframe format 

str(Uptake.ratio.np.all)


str(Uptake.ratio.gg)

Uptake.ratio.gg.all<- Uptake.ratio.gg[,1:8] %>% tidyr::gather(sample, ratio, 3:8) #transform from wide to long dataframe format 



########################################    plot uptake ratios maps whole genome      #######################################



u<- unique(Uptake.ratio.np.all$sample)

for (i in 1:length(u)){
  p<- Uptake.ratio.np.all %>% dplyr::filter(sample == u[i]) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 20), expand = c(0, 0))+
    labs(x = "NP genome positions", y = "uptake ratio") +
    ggtitle(" coverage maps NP donor DNA") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/uptake_maps","whole_genome", u[i],"2", "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}


str(Uptake.ratio.gg.all)

u<- unique(Uptake.ratio.gg.all$sample)

for (i in 4:length(u)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u[i]) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 30), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "uptake ratio") +
    ggtitle(" coverage maps PittGG donor DNA") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/uptake_maps","whole_genome","2", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}




########################################    plot uptake ratios maps 1kb, 10kb, 100kb long donor DNA data      #######################################



u<- unique(Uptake.ratio.np.all$sample)

for (i in 1:length(u)){
  p<- Uptake.ratio.np.all %>% dplyr::filter(sample == u[i]) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1000, 100), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 20), expand = c(0, 0))+
    labs(x = "NP genome positions", y = "uptake ratio") +
    ggtitle(" coverage maps NP donor DNA") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/long_1kb/uptake_maps","1kb", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}


str(Uptake.ratio.gg.all)

u<- unique(Uptake.ratio.gg.all$sample)

for (i in 4:length(u)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u[i]) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 30), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "uptake ratio") +
    ggtitle(" coverage maps PittGG donor DNA") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/uptake_maps","whole_genome","2", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}