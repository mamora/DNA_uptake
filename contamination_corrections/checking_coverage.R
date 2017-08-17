

#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
#########################################################################

#############################load dataframes##############################
depth.rdnp.map<- fread("~/DNA_uptake/datasets/tables/rdnp.map.csv") #load concataneted samples aligned to concatenated rdnp genomes and only including reads with quality > 0  
depth.rdnp.all<- fread("~/DNA_uptake/datasets/tables/rdnp.all.f.csv") #load concataneted samples aligned to concatenated rdnp genomes including all reads 
colnames(depth.rdnp.map)<- c("name","genome","pos","depth") # add headers
colnames(depth.rdnp.all)<- c("name","genome","pos","depth") # add headers
#########################################################################

#############################check dataframes##############################
str(depth.rdnp.map)
str(depth.rdnp.all)
names<- unique(depth.rdnp.map$name) #check sample names
unique(depth.rdnp.all$name) #check sample names
#########################################################################


#############################save or load dataframes (if needed)##############################
write.csv(UP07.depth.rdnp, "~/DNA_uptake/datasets/tables/UP07.depth.rdnp.csv") #
##############################################################################################





##############################plot coverage rdnp all first 1kb all samples KW20 #####################################

theme_set(theme_grey()) #set theme

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_all/KW20/" # Where these files are located
# Set directory
setwd(whereami)

plot.coverage.(dat = depth.rdnp.all, t = "all",gen = "KW20", size = 10000) 

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_all/KW20_2/" # Where these files are located
# Set directory
setwd(whereami)

plot.coverage.part(dat = depth.rdnp.all, t = "all",gen = "KW20", start = 193000, end = 195000) 



### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_all/KW20all/" # Where these files are located
# Set directory
setwd(whereami)

plot.coverage.all(dat = depth.rdnp.all, t = "all",gen = "KW20", size = 1831585)
                 
                   


for(i in 1:length(names)){
  rd1<- dplyr::filter(depth.rdnp.all, name == names[i])
  u<- rd1 %>%  dplyr::filter(genome == "KW20") %>% slice(1:10000) %>% 
    ggplot(aes(x = pos,y = depth)) +
    geom_point(shape = 20, size = 1)+
    scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 600))+
    labs(x = "KW20 genome positions", y = "depth") +
    ggtitle(" raw coverage rdnp aligned to KW20, all reads") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("coverage_depth", names[i], ".tiff", sep="_")
  tiff(file_name, width = 800, height = 500, units = "px")
  print(u)
  dev.off()  
}

#############################################################################################################

  

  
  
##############################plot coverage rdnp all first 1kb all samples NP #####################################

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_all/NP/" # Where these files are located

# Set directory
setwd(whereami)

plot.coverage.(dat = depth.rdnp.all, t = "all",gen = "NP", size = 10000) 

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_all/NP_2/" # Where these files are located
# Set directory
setwd(whereami)

plot.coverage.part(dat = depth.rdnp.all, t = "all",gen = "NP", start = 192000, end = 196000) 





### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_all/NPall/" # Where these files are located
# Set directory
setwd(whereami)

plot.coverage.all(dat = depth.rdnp.all, t = "all",gen = "NP", size = 1914386) 



for(i in 1:length(names)){
  rd1<- dplyr::filter(depth.rdnp.all, name == names[i])
  u<- rd1 %>%  dplyr::filter(genome == "NP") %>% slice(1:10000) %>% 
    ggplot(aes(x = pos,y = depth)) +
    geom_point(shape = 20, size = 1)+
    scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 600))+
    labs(x = "NP genome positions", y = "depth") +
    ggtitle(" raw coverage rdnp aligned to NP, all reads") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("coverage_depth", names[i], ".tiff", sep="_")
  tiff(file_name, width = 800, height = 500, units = "px")
  print(u)
  dev.off()  
}

#############################################################################################################






##############################plot coverage rdnp all first 1kb map samples KW20 #####################################

theme_set(theme_grey()) #set theme

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_map/KW20/" # Where these files are located

# Set directory
setwd(whereami)

plot.coverage.(dat = depth.rdnp.map, t = "Map > 0",gen = "KW20", size = 10000) 

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_map/KW20all/" # Where these files are located

# Set directory
setwd(whereami)

plot.coverage.all(dat = depth.rdnp.map, t = "Map > 0",gen = "KW20", size = 1831585) 




for(i in 1:length(names)){
  rd1<- dplyr::filter(depth.rdnp.map, name == names[i])
  u<- rd1 %>%  dplyr::filter(genome == "KW20") %>% slice(1:10000) %>% 
    ggplot(aes(x = pos,y = depth)) +
    geom_point(shape = 20, size = 1)+
    scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 600))+
    labs(x = "KW20 genome positions", y = "depth") +
    ggtitle(" raw coverage rdnp aligned to KW20, map reads") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("coverage_depth", names[i], ".tiff", sep="_")
  tiff(file_name, width = 800, height = 500, units = "px")
  print(u)
  dev.off()  
}

#############################################################################################################




 
##############################plot coverage rdnp all first 1kb all samples NP #####################################

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_map/NP/" # Where these files are located

# Set directory
setwd(whereami)


plot.coverage.(dat = depth.rdnp.map, t = "Map > 0",gen = "NP", size = 10000) 

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_map/NP_2/" # Where these files are located
# Set directory
setwd(whereami)

plot.coverage.part(dat = depth.rdnp.map, t = "Map > 0",gen = "NP", start = 192000, end = 196000) 



### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/task_list_marcelo/correccion/coveragegraphs/rdnp_map/NPall/" # Where these files are located

# Set directory
setwd(whereami)

plot.coverage.all(dat = depth.rdnp.map, t = "Map > 0",gen = "NP", size = 1914386) 




plot.coverage.part<- function(dat = depth.rdnp.map, t = "Map > 0",gen = "NP", start = 1, end = 10000){
for(i in 1:length(names)){
  rd1<- dplyr::filter(dat, name == names[i])
  u<- rd1 %>%  dplyr::filter(genome == gen) %>% slice(start:end) %>% 
    ggplot(aes(x = pos,y = depth)) +
    geom_point(shape = 20, size = 1)+
    scale_x_continuous(breaks = seq(start, end,  500), expand = c(0, 0))+
    scale_y_continuous(limits = c(0, 600))+
    labs(x = paste(gen,"genome positions", sep=''), y = "depth") +
    ggtitle(paste(names[i], 'raw coverage \n', 
                  "aligned to", "NP",
                  t,
                  sep=' ')) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  file_name = paste("coverage_depth", names[i], ".tiff", sep="_")
  tiff(file_name, width = 800, height = 500, units = "px")
  print(u)
  dev.off()  
}
}

plot.coverage.all<- function(dat = depth.rdnp.map, t = "Map > 0",gen = "NP", size = 10000){
  for(i in 1:length(names)){
    rd1<- dplyr::filter(dat, name == names[i])
    u<- rd1 %>%  dplyr::filter(genome == gen) %>% slice(1:size) %>% 
      ggplot(aes(x = pos,y = depth)) +
      geom_point(shape = ".", size = 1)+
      scale_x_continuous(breaks = seq(0 , size, 200000), expand = c(0, 0))+
      scale_y_continuous(limits = c(0, 600))+
      labs(x = paste(gen,"genome positions", sep=''), y = "depth") +
      ggtitle(paste(names[i], 'raw coverage \n', 
                    "aligned to", "NP",
                    t,
                    sep=' ')) +
      theme(plot.margin=unit(c(1,1,1,1),"cm"),
            legend.position = "bottom",
            panel.grid.minor = element_line(colour="white", size=0.5),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.text  = element_text(size=18),
            axis.title = element_text(size = 18, face = "bold")) 
    file_name = paste("coverage_depth", names[i], ".tiff", sep="_")
    tiff(file_name, width = 1000, height = 500, units = "px")
    print(u)
    dev.off()  
  }
}

#############################################################################################################










