
#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
#########################################################################


#############################load dataframes##############################
input.all<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/input.all.csv") #load samples
rd<-  fread("~/DNA_uptake/datasets/RR722.rd.pilon.all.depth.bed") 
str(rd)
unique(rd$V4)
rd.all<- rd %>% filter(V4 == "gi|16271976|ref|NC_000907.1|pilon")
colnames(input.all)<- c("name","genome","pos","depth") # add headers
u13<- dplyr::filter(input.all, name == "UP13")
u14<- dplyr::filter(input.all, name == "UP14")
u15<- dplyr::filter(input.all, name == "UP15")
u16<- dplyr::filter(input.all, name == "UP16")

#########################################################################

#############################check dataframes##############################
str(input.all)
u<- unique(input.all$name) #check sample names
#########################################################################


#1 kb, 10 kb, 100 kb and whole-genome maps for UP13-16 without added Rd.  Comparing these will let us distinguish between random variation and sequencing biases. Use normalized coverages for these graphs, for easier comparisons


##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up13.all.mapped<- 2723796 # total mapped reads from summary table
up13.n<- (u13$depth *1e+6)/up13.all.mapped
up14.all.mapped<- 2812556 # total mapped reads from summary table
up14.n<- (u14$depth *1e+6)/up14.all.mapped
up15.all.mapped<- 4750647 # total mapped reads from summary table
up15.n<- (u15$depth *1e+6)/up15.all.mapped
up16.all.mapped<- 10168072 # total mapped reads from summary table
up16.n<- (u16$depth *1e+6)/up16.all.mapped
#########################################################################


#############################   make normalized coverage dataframes   #######################

input.n.np<- data.frame(pos = u13$pos, up13.n = up13.n, up15.n =up15.n) 

input.n.gg<- data.frame(pos = u14$pos, up14.n = up14.n,  up16.n = up16.n) 
#############################################################################################

#######################    tidy up dataframes   ##########################
input.n.np.t<- input.n.np %>% tidyr::gather(sample, "n.depth", 2:3)
unique(input.n.np.t$sample)
input.n.gg.t<- input.n.gg %>% tidyr::gather(sample, "n.depth", 2:3)
unique(input.n.gg.t$sample)

##########################################################################



##################################    Plot normalized coverage maps #####################################


u<- input.n.np %>%  dplyr::slice(1:1000) %>% tidyr::gather(sample, "n.depth", 2:3) %>% 
  ggplot(aes(x = pos,y = n.depth)) +
  geom_point(aes(color = sample), shape = 20, size = 1)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_x_continuous(breaks = seq(0 , 1000, 100), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 800))+
  labs(x = "NP genome positions", y = "normalized depth") +
  ggtitle(" coverage map normalized input ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "1kb","NP", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off() 


u<- input.n.np %>%  dplyr::slice(1:10000) %>% tidyr::gather(sample, "n.depth", 2:3) %>% 
  ggplot(aes(x = pos,y = n.depth)) +
  geom_point(aes(color = sample), shape = ".", size = 1)+
  guides(colour = guide_legend(override.aes = list(shape=20, size = 3))) +  
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 800))+
  labs(x = "NP genome positions", y = "normalized depth") +
  ggtitle(" coverage map normalized input ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "10kb","2","NP", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off()  



u<- input.n.np %>%  dplyr::slice(1:100000) %>% tidyr::gather(sample, "n.depth", 2:3) %>% 
  ggplot(aes(x = pos,y = n.depth)) +
  geom_point(aes(color = sample), shape = ".", size = 1)+
  guides(colour = guide_legend(override.aes = list(shape=20, size = 3))) +
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 800))+
  labs(x = "NP genome positions", y = "normalized depth") +
  ggtitle(" coverage map normalized input ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "100kb","NP", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off()  



u<- input.n.gg %>%  dplyr::slice(1:100000) %>% tidyr::gather(sample, "n.depth", 2:3) %>% 
  ggplot(aes(x = pos,y = n.depth)) +
  geom_point(aes(color = sample), shape = ".", size = 1)+
  guides(colour = guide_legend(override.aes = list(shape=20, size = 3))) +
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 800))+
  labs(x = "PittGG genome positions", y = "normalized depth") +
  ggtitle(" coverage map normalized input ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "100kb","PittGG", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off() 




u<- input.n.gg %>%  dplyr::slice(1:10000) %>% tidyr::gather(sample, "n.depth", 2:3) %>% 
  ggplot(aes(x = pos,y = n.depth)) +
  geom_point(aes(color = sample), shape = ".", size = 1)+
  guides(colour = guide_legend(override.aes = list(shape=20, size = 3))) +
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 800))+
  labs(x = "PittGG genome positions", y = "normalized depth") +
  ggtitle(" coverage map normalized input ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "10kb","PittGG", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off() 



u<- input.n.gg %>%  dplyr::slice(10000:11000) %>% tidyr::gather(sample, "n.depth", 2:3) %>% 
  ggplot(aes(x = pos,y = n.depth)) +
  geom_point(aes(color = sample), shape = 20, size = 1)+
  guides(colour = guide_legend(override.aes = list(shape=20, size = 3))) +
  scale_x_continuous(breaks = seq(10000 , 11000, 100), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 800))+
  labs(x = "PittGG genome positions", y = "normalized depth") +
  ggtitle(" coverage map normalized input ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "1kb","PittGG", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off() 


######################################    Plot histogram   #####################################




p <- ggplot() +
  geom_histogram(aes(x = depth), binwidth = 25, colour = "black", data = input.all) +
  facet_wrap(~name) +
  labs(x = "coverage") +
  ggtitle("Histogram of coverage of input samples without correction") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/coverage_hist","no_correction", "tiff", sep=".")
tiff(file_name, width = 1600, height = 900, units = "px")
print(p)
dev.off()

