
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
rd.all.mapped<- 4088504 # total mapped reads from summary table
rd.n<- (rd.all$depth *1e+6)/rd.all.mapped
#########################################################################

#############################   make normalized coverage dataframes   #######################

rd.n.all<- data.frame(pos = rd.all$pos, depth.n = rd.n) 

write.csv(rd.n.all, "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/datasets/rd.n.all.csv")

#############################################################################################


##################################    Plot normalized coverage maps #####################################

#  1.  Positions 1187001-1191000


u<- rd.n.all %>%  dplyr::slice(1187001:1191000) %>% 
  ggplot(aes(x = pos,y = depth.n)) +
  geom_point(shape = ".", size = 1)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 500), expand = c(0, 0))+
  labs(x = "Rd genome positions", y = "normalized depth") +
  ggtitle(" Rd coverage map normalized depth ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "1.1","Rd", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off() 


#  2.  Positions 1734501-1738000


u<- rd.n.all %>%  dplyr::slice(1734501:1738000) %>% 
  ggplot(aes(x = pos,y = depth.n)) +
  geom_point(shape = ".", size = 1)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 500), expand = c(0, 0))+
  labs(x = "Rd genome positions", y = "normalized depth") +
  ggtitle(" Rd coverage map normalized depth ") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/coverage/norm_input_map/normalized_input", "2","Rd", "tiff", sep=".")
tiff(file_name, width = 900, height = 500, units = "px")
print(u)
dev.off()