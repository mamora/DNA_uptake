
#####################################    Goal     ######################################### 

#To plot coverage of input and uptake samples as well as coverage histograms. 


#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
#########################################################################

#############################load dataframes##############################
##WARNING MEMORY INTENSIVE###  
#YOU MIGHT HAVE TO LOAD ONE AT A TIME AND REFRESH THE SESSION OFTEN
# load dataframes
raw.depth.samples<- fread("~/DNA_uptake/datasets/tables/raw.depth.samples.csv") #load new input samples fro UP01   
raw.depth.inputs<- fread("~/DNA_uptake/datasets/tables/raw.depth.inputs.csv") #load new input samples fro UP01   
Uptake.ratio.np<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.np.csv") #load new input samples fro UP01   
# load dataframes
raw.depth.samples.gg<- fread("~/DNA_uptake/datasets/tables/raw.depth.samples.gg.csv") #load new input samples fro UP01   
raw.depth.inputs.gg<- fread("~/DNA_uptake/datasets/tables/raw.depth.inputs.gg.csv") #load new input samples fro UP01   
Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv") #load new input samples fro UP01   

#########################################################################



###############################    turn wide dataframes into long dataframes #########################
######################################################################################################


theme_set(theme_grey()) #set theme

str(raw.depth.inputs) #check dataframe

str(raw.depth.samples)


# let's tidy up things a bit!
raw.depth.inputs.all<- raw.depth.inputs %>% tidyr::gather(sample, coverage, 3:10) #transform from wide to long dataframe format 

raw.depth.samples.all<- raw.depth.samples %>% tidyr::gather(sample, coverage, 3:8) #transform from wide to long dataframe format 

raw.depth.inputs.gg.all<- raw.depth.inputs.gg %>% tidyr::gather(sample, coverage, 3:10) #transform from wide to long dataframe format 

raw.depth.samples.gg.all<- raw.depth.samples.gg %>% tidyr::gather(sample, coverage, 3:8) #transform from wide to long dataframe format 



str(raw.depth.inputs.all) #check dataframe

str(raw.depth.samples.all)

str(raw.depth.inputs.gg.all) #check dataframe

str(raw.depth.samples.gg.all)

unique(raw.depth.inputs.all$sample) #check samples


###############################  Plot Coverage histograms ##############################


p <- ggplot() +
  geom_histogram(aes(x = coverage), binwidth = 25, colour = "black", data = raw.depth.inputs.all) +
  labs(x = "coverage") +
  facet_wrap(~sample) +
  ggtitle("Histogram of coverage of input samples for NP uptake samples") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_hist","input","NP", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = coverage), binwidth = 25, colour = "black", data = raw.depth.samples.all) +
  labs(x = "coverage") +
  facet_wrap(~sample) +
  ggtitle("Histogram of coverage of samples samples for NP uptake samples") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_hist","uptake","NP", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()





p <- ggplot() +
  geom_histogram(aes(x = coverage), binwidth = 25, colour = "black", data = raw.depth.inputs.gg.all) +
  labs(x = "coverage") +
  facet_wrap(~sample) +
  ggtitle("Histogram of coverage of input samples for PittGG uptake samples") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_hist","input","PittGG", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


p <- ggplot() +
  geom_histogram(aes(x = coverage), binwidth = 25, colour = "black", data = raw.depth.samples.gg.all) +
  labs(x = "coverage") +
  facet_wrap(~sample) +
  ggtitle("Histogram of coverage of samples samples for PittGG uptake samples") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_hist","uptake","PittGG", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()


u<- unique(raw.depth.samples.all$sample)

for (i in 1:length(u)){
  p<- raw.depth.samples.all %>% dplyr::filter(sample == u[i]) %>% ggplot() +
      geom_histogram(aes(x = coverage), binwidth = 25, colour = "black") +
     scale_x_continuous(limits = c(0,1000), expand = c(0, 0))+  
     scale_y_continuous(limits = c(0,7e+5), expand = c(0, 0))+  
     labs(x = "coverage") +
      ggtitle(paste("Histogram of coverage from ", u[i], "NP","uptake donor DNA", sep=" ")) +
      theme(plot.margin=unit(c(1,1,1,1),"cm"),
            legend.position = "bottom",
            panel.grid.minor = element_line(colour="white", size=0.5),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.text  = element_text(size=18),
            axis.title = element_text(size = 18, face = "bold")) 
    file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_hist","uptake", u[i], "tiff", sep=".")
    tiff(file_name, width = 1000, height = 500, units = "px")
    print(p)
    dev.off()
}  



u<- unique(raw.depth.samples.gg.all$sample)

for (i in 1:length(u)){
  p<- raw.depth.samples.gg.all %>% dplyr::filter(sample == u[i]) %>% ggplot() +
    geom_histogram(aes(x = coverage), binwidth = 25, colour = "black") +
    scale_x_continuous(limits = c(0,1000), expand = c(0, 0))+  
    scale_y_continuous(limits = c(0,7e+5), expand = c(0, 0))+  
    labs(x = "coverage") +
    ggtitle(paste("Histogram of coverage from ", u[i], "GG","uptake donor DNA", sep=" ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_hist","uptake", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 500, units = "px")
  print(p)
  dev.off()
}



#########################################    Plot input coverage maps    #######################3

str(raw.depth.inputs.all)

u<- unique(raw.depth.inputs.all$sample)

for (i in 1:length(u)){
  p<- raw.depth.inputs.all %>% dplyr::filter(sample == u[i]) %>% 
  ggplot(aes(x = pos,y = coverage)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 600), expand = c(0, 0))+
    labs(x = "NP genome positions", y = "Depth") +
    ggtitle(" coverage maps NP donor DNA") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_maps","whole_genome", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}



str(raw.depth.inputs.gg.all)

u<- unique(raw.depth.inputs.gg.all$sample)

for (i in 1:length(u)){
  p<- raw.depth.inputs.gg.all %>% dplyr::filter(sample == u[i]) %>% 
    ggplot(aes(x = pos,y = coverage)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 600), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "Depth") +
    ggtitle(" coverage maps PittGG donor DNA") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_maps","whole_genome", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}




#################################   plot coverage maps of uptake samples ##############################


#########################################    Plot uptake samples coverage maps    #######################3

str(raw.depth.samples.all)

u<- unique(raw.depth.samples.all$sample)

for (i in 1:length(u)){
  p<- raw.depth.samples.all %>% dplyr::filter(sample == u[i]) %>% 
    ggplot(aes(x = pos,y = coverage)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 5000), expand = c(0, 0))+
    labs(x = "NP genome positions", y = "Depth") +
    ggtitle(paste("coverage maps", u[i], "donor DNA", sep=".")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_maps","whole_genome", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}

str(raw.depth.samples.all)

u<- unique(raw.depth.samples.gg.all$sample)

for (i in 1:length(u)){
  p<- raw.depth.samples.gg.all %>% dplyr::filter(sample == u[i]) %>% 
    ggplot(aes(x = pos,y = coverage)) +
    geom_point(shape = ".", size = 1)+
    scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 2000), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "Depth") +
    ggtitle(paste("coverage maps", u[i], "donor DNA", sep=".")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_maps_uptake","whole_genome", u[i],"2", "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(p)
  dev.off()  
}



#################################   Plot uptake ratios histograms     #####################################



str(Uptake.ratio.np.all)


for (i in 1:length(u)){
  p<- Uptake.ratio.np.all %>% dplyr::filter(sample == u[i]) %>% ggplot() +
    geom_histogram(aes(x = ratio), binwidth = 0.5, colour = "black") +
    scale_y_continuous(limits = c(0,1.2e+6), expand = c(0, 0))+  
    scale_x_continuous(limits = c(0,20), breaks = seq(0 , 20, 2), expand = c(0, 0))+ 
    labs(x = "uptake ratios") +
    ggtitle(paste("Histogram of uptake ratios from ", u[i],"uptake donor DNA", sep=" ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/uptake_ratios_histograms/histogram","ratios", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 500, units = "px")
  print(p)
  dev.off()
}  



for (i in 1:length(u)){
  p<- Uptake.ratio.gg.all %>% dplyr::filter(sample == u[i]) %>% ggplot() +
    geom_histogram(aes(x = ratio), binwidth = 0.5, colour = "black") +
    scale_y_continuous(limits = c(0,1.2e+6), expand = c(0, 0))+  
    scale_x_continuous(limits = c(0,20), breaks = seq(0 , 20, 2), expand = c(0, 0))+ 
    labs(x = "uptake ratios") +
    ggtitle(paste("Histogram of uptake ratios from ", u[i],"uptake donor DNA", sep=" ")) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/uptake_ratios_histograms/histogram","ratios", u[i], "tiff", sep=".")
  tiff(file_name, width = 1000, height = 500, units = "px")
  print(p)
  dev.off()
}
