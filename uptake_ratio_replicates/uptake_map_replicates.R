
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
  p<- Uptake.ratio.np.all %>% dplyr::filter(sample == u[i]) %>% dplyr::slice(1:1000) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1000, 100), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
    labs(x = "NP genome positions", y = "uptake ratio") +
    ggtitle(paste("uptake ratio map", u[i], "NP donor DNA", sep = " ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/1kb/uptake_maps","1kb", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}


for (i in 1:length(u)){
  p<- Uptake.ratio.np.all %>% dplyr::filter(sample == u[i]) %>% dplyr::slice(1:10000) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
    labs(x = "NP genome positions", y = "uptake ratio") +
    ggtitle(paste("uptake ratio map", u[i], "NP donor DNA", sep = " ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/10kb/uptake_maps","10kb", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}

for (i in 1:length(u)){
  p<- Uptake.ratio.np.all %>% dplyr::filter(sample == u[i]) %>% dplyr::slice(1:100000) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
    labs(x = "NP genome positions", y = "uptake ratio") +
    ggtitle(paste("uptake ratio map", u[i], "NP donor DNA", sep = " ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/100kb/uptake_maps","100kb", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}

u1<- unique(Uptake.ratio.gg.all$sample)

for (i in 1:length(u1)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u1[i]) %>% dplyr::slice(1:1000) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1000, 100), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 25), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "uptake ratio") +
    ggtitle(paste("uptake ratio map", u1[i], "PittGG donor DNA", sep = " ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/1kb/uptake_maps","1kb", u1[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}



for (i in 1:length(u1)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u1[i]) %>% dplyr::slice(1:10000) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 25), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "uptake ratio") +
    ggtitle(paste("uptake ratio map", u1[i], "PittGG donor DNA", sep = " ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/10kb/uptake_maps","10kb", u1[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}


for (i in 1:length(u1)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u1[i]) %>% dplyr::slice(1:100000) %>% 
    ggplot(aes(x = pos,y = ratio)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 25), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "uptake ratio") +
    ggtitle(paste("uptake ratio map", u1[i], "PittGG donor DNA", sep = " ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/replicates_uptake/100kb/uptake_maps","100kb", u1[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}


########################################################################################################################
#################################   Run correlations between long and short uptake ratios   ############################

str(Uptake.ratio.np)


corr<- cor(Uptake.ratio.np$ratio_short, Uptake.ratio.np$ratio_long)



p<- Uptake.ratio.np %>% dplyr::slice(1:100000) %>% 
  ggplot(aes(x = ratio_short,y = ratio_long)) +
  geom_point(shape = ".", size = 1)+
  geom_text(x = 2 , y = 3 ,label = "correlation = 0.28" ) +
  scale_x_continuous(limits = c(0, 6), expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "NP short uptake ratio", y = "NP long uptake ratio") +
  ggtitle("uptake ratio of short vs long") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/uptake_ratio_corr_NP","100k", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(p)
dev.off()  

p<- Uptake.ratio.np  %>% 
  ggplot(aes(x = ratio_short,y = ratio_long)) +
  geom_point(shape = ".", size = 1)+
  geom_text(x = 2 , y = 3 ,label = "correlation = 0.28" ) +
  scale_x_continuous(limits = c(0, 6), expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "NP short uptake ratio", y = "NP long uptake ratio") +
  ggtitle("uptake ratio of short vs long") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/uptake_ratio_corr_NP","whole", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(p)
dev.off()


################   remove infinite and missing values to calculate the correlation ###############3 
a<- which(Uptake.ratio.gg$ratio_long == Inf)

a1<- which(is.na(Uptake.ratio.gg$ratio_long) == TRUE)

a2<- which(is.na(Uptake.ratio.gg$ratio_short) == TRUE)

u<- unique(c(a,a1,a2)) # get all the numbers without repetitions


s<- Uptake.ratio.gg$ratio_short[-u]  #remove missing and infinite

s1<- Uptake.ratio.gg$ratio_long[-u]

max(s1)  #what??? somewhere with uptake ratio of 175??

cor(s, s1)

str(Uptake.ratio.gg)
 

p<- Uptake.ratio.gg  %>% 
  ggplot(aes(x = ratio_short,y = ratio_long)) +
  geom_point(shape = ".", size = 1)+
  geom_text(x = 8 , y = 8.5 ,label = "correlation = 0.24" ) +
  scale_x_continuous(limits = c(0, 30), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG short uptake ratio", y = "PittGG long uptake ratio") +
  ggtitle("uptake ratio of short vs long") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/uptake_ratio_corr_GG","whole", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(p)
dev.off()


