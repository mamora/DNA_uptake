
#####################################    Goal     ######################################### 

#To generate uptake ratios for small and large fragment data for PittGG donor DNA samples




#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
#########################################################################

#############################load dataframes##############################
remix14.UP04<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX14.UP04.I.gg.depth.bed") #load new input samples fro UP01   
colnames(remix14.UP04)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
samples.gg<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/samples.all.gg.bed.csv") #load all samples from NP donor DNA, input included
colnames(samples.gg)<- c("sample","genome","pos","depth") # add headers
up04<- dplyr::filter(samples.gg, sample == "UP04") #subset UP01 
remix14.UP04$UP04<- up04$depth
remix14.UP04$input.mean<- (remix14.UP04$depth.1+remix14.UP04$depth.2+remix14.UP04$depth.3+remix14.UP04$depth.4)/4
#########################################################################

n<- length(which(remix14.UP04$input.mean <= 10))

(n/length(remix14.UP04$input.mean)*100)


#############################check dataframes##############################
str(samples.gg)
str(remix14.UP04)
unique(samples.gg$sample) #check sample names
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization#### ##########################
up14.up04.mapped<-  sum(remix14.UP04$input.mean)/length(remix14.UP04$input.mean)  
up04.all.mapped<- sum(remix14.UP04$UP04)/length(remix14.UP04$input.mean)
remix14.UP04$UP04.n<- (remix14.UP04$UP04)/up04.all.mapped
remix14.UP04$input.n<- (remix14.UP04$input.mean)/up14.up04.mapped
remix14.UP04$UP04.ratio<- remix14.UP04$UP04.n/remix14.UP04$input.n
#########################################################################

g<- Uptake.ratio.gg$ratio_short[-(which(Uptake.ratio.gg$ratio_short == Inf))]

mean(g, na.rm = TRUE)

h<- mean(remix14.UP04$UP04.n)

##############################plot coverage rdnp all first 1kb all samples KW20 #####################################

theme_set(theme_grey()) #set theme

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/samples_figs/"
# Set directory
setwd(whereami)


u<- remix14.UP04 %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = UP04.ratio)) +
  geom_point(shape = 20, size = 1)+
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP04") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("UP04.ratio","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  


u<- remix14.UP04 %>% 
  ggplot(aes(x = pos,y = UP04.ratio)) +
  geom_point(shape = ".", size = 1)+
  scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP04") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))
file_name = paste("UP04.ratio","whole_genome", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off() 


###################################     UP05       ######################################

#############################load dataframes##############################
remix14.UP05<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX14.UP05.I.gg.depth.bed") #load new input samples for UP05   
colnames(remix14.UP05)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up05<- dplyr::filter(samples.gg, sample == "UP05") #subset UP05 
remix14.UP05$UP05<- up05$depth
remix14.UP05$input.mean<- (remix14.UP05$depth.1+remix14.UP05$depth.2+remix14.UP05$depth.3+remix14.UP05$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(samples.np)
str(remix14.UP05)
#########################################################################



##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up14.up05.mapped<- sum(remix14.UP05$input.mean)/length(remix14.UP04$input.mean)   # total mapped reads, from sambamba flagstats files 
up05.all.mapped<- sum(remix14.UP05$UP05)/length(remix14.UP04$input.mean)   # total mapped reads from summary table
remix14.UP05$UP05.n<- (remix14.UP05$UP05)/up05.all.mapped
remix14.UP05$input.n<- (remix14.UP05$input.mean)/up14.up05.mapped
remix14.UP05$UP05.ratio<- remix14.UP05$UP05.n/remix14.UP05$input.n
#########################################################################

mean(remix14.UP05$UP05.n)

sum(remix14.UP05$UP05.ratio)/length(remix14.UP05$UP05.ratio)

which(remix14.UP05$UP05.ratio == Inf)

g<- remix14.UP05$UP05.ratio[-which(remix14.UP05$UP05.ratio == Inf)]

mean(g, na.rm = TRUE)

#################################################################################################

###################################     UP06       ######################################

#############################load dataframes##############################
remix14.UP06<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX14.UP06.I.gg.depth.bed") #load new input samples for UP06   
colnames(remix14.UP06)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up06<- dplyr::filter(samples.gg, sample == "UP06") #subset UP06 
remix14.UP06$UP06<- up06$depth
remix14.UP06$input.mean<- (remix14.UP06$depth.1+remix14.UP06$depth.2+remix14.UP06$depth.3+remix14.UP06$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(remix14.UP06)
#########################################################################



##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up14.up06.mapped<- sum(remix14.UP06$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up06.all.mapped<- sum(remix14.UP06$UP06)/length(remix14.UP06$input.mean)  # total mapped reads from summary table
remix14.UP06$UP06.n<- (remix14.UP06$UP06)/up06.all.mapped
remix14.UP06$input.n<- (remix14.UP06$input.mean)/up14.up06.mapped
remix14.UP06$UP06.ratio<- remix14.UP06$UP06.n/remix14.UP06$input.n
#########################################################################

####################################     save in dataframes     ###############################
raw.depth.samples.gg<- data.table(pos = remix14.UP04$pos, UP04 = remix14.UP04$UP04, UP05 = remix14.UP05$UP05, UP06 = remix14.UP06$UP06) 
raw.depth.inputs.gg<-  data.table(pos = remix14.UP04$pos, UP14.UP04 = remix14.UP04$input.mean, UP14.UP05 = remix14.UP05$input.mean, UP14.UP06 = remix14.UP06$input.mean) 
Uptake.ratio.gg<-  data.table(pos = remix14.UP04$pos, UP04.ratio = remix14.UP04$UP04.ratio, UP05.ratio = remix14.UP05$UP05.ratio, UP06.ratio = remix14.UP06$UP06.ratio) 




#################################      UP07       ############################################

#############################load dataframes##############################
remix16.UP10<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX16.UP10.I.gg.depth.bed") #load new input samples fro UP10   
colnames(remix16.UP10)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up10<- dplyr::filter(samples.gg, sample =="UP10") #subset UP10 
remix16.UP10$UP10<- up10$depth
remix16.UP10$input.mean<- (remix16.UP10$depth.1+remix16.UP10$depth.2+remix16.UP10$depth.3+remix16.UP10$depth.4)/4

#########################################################################


#############################check dataframes##############################
str(remix16.UP10)
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up16.up10.mapped<- sum(remix16.UP10$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up10.all.mapped<- sum(remix16.UP10$UP10)/length(remix14.UP06$input.mean) # total mapped reads from summary table
remix16.UP10$UP10.n<- (remix16.UP10$UP10)/up10.all.mapped
remix16.UP10$input.n<- (remix16.UP10$input.mean)/up16.up10.mapped
remix16.UP10$UP10.ratio<- remix16.UP10$UP10.n/remix16.UP10$input.n
#########################################################################

sum(remix16.UP10$UP10)
sum(remix16.UP10$input.mean)

# save in dataframes
raw.depth.samples.gg$UP10<- remix16.UP10$UP10
raw.depth.inputs.gg$UP16.UP10<- remix16.UP10$input.mean
Uptake.ratio.gg$UP10.ratio<- remix16.UP10$UP10.ratio


##############################plot coverage rdnp all first 1kb all samples KW20 #####################################

theme_set(theme_grey()) #set theme

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/samples_figs/"
# Set directory
setwd(whereami)


u<- remix16.UP10 %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = UP10.ratio)) +
  geom_point(shape = 20, size = 1)+
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP10") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("UP10.ratio.2","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off() 


u<- remix16.UP10 %>% 
  ggplot(aes(x = pos,y = UP10.ratio)) +
  geom_point(shape = ".", size = 1)+
  scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP10") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))
file_name = paste("UP10.ratio","whole_genome", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off() 


############################################################################################

######################################      UP11        ###################################################



#############################load dataframes##############################
remix16.UP11<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX16.UP11.I.gg.depth.bed") #load new input samples for UP11   
colnames(remix16.UP11)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up11<- dplyr::filter(samples.gg, sample =="UP11") #subset UP07 
remix16.UP11$UP11<- up11$depth
remix16.UP11$input.mean<- (remix16.UP11$depth.1+remix16.UP11$depth.2+remix16.UP11$depth.3+remix16.UP11$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(remix16.UP11)
str(samples.gg)
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up16.up11.mapped<- sum(remix16.UP11$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up11.all.mapped<- sum(remix16.UP11$UP11)/length(remix14.UP06$input.mean)  # total mapped reads from summary table
remix16.UP11$UP11.n<- (remix16.UP11$UP11)/up11.all.mapped
remix16.UP11$input.n<- (remix16.UP11$input.mean)/up16.up11.mapped
remix16.UP11$UP11.ratio<- remix16.UP11$UP11.n/remix16.UP11$input.n
#########################################################################

sum(remix16.UP11$UP11)
sum(remix16.UP11$input.mean)


raw.depth.samples.gg$UP11<- remix16.UP11$UP11
raw.depth.inputs.gg$UP16.UP11<- remix16.UP11$input.mean
Uptake.ratio.gg$UP11.ratio<- remix16.UP11$UP11.ratio

#######################################################################################################


######################################      UP12        ###################################################



#############################load dataframes##############################
remix16.UP12<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX16.UP12.I.gg.depth.bed") #load new input samples fro UP01   
colnames(remix16.UP12)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up12<- dplyr::filter(samples.gg, sample =="UP12") #subset UP12 
remix16.UP12$UP12<- up12$depth
remix16.UP12$input.mean<- (remix16.UP12$depth.1+remix16.UP12$depth.2+remix16.UP12$depth.3+remix16.UP12$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(remix16.UP12)
str(samples.gg)
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up16.up12.mapped<- sum(remix16.UP12$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up12.all.mapped<- sum(remix16.UP12$UP12)/length(remix14.UP06$input.mean) # total mapped reads from summary table
remix16.UP12$UP12.n<- (remix16.UP12$UP12)/up12.all.mapped
remix16.UP12$input.n<- (remix16.UP12$input.mean)/up16.up12.mapped
remix16.UP12$UP12.ratio<- remix16.UP12$UP12.n/remix16.UP12$input.n
#########################################################################

raw.depth.samples.gg$UP12<- remix16.UP12$UP12
raw.depth.inputs.gg$UP16.UP12<- remix16.UP12$input.mean
Uptake.ratio.gg$UP12.ratio<- remix16.UP12$UP12.ratio

#######################################################################################################



##################calculate mean and standard deviation of the three replicates#######################

Uptake.ratio.gg$ratio_long<- apply(Uptake.ratio.gg[,2:4], 1, mean)

Uptake.ratio.gg$ratio_short<- apply(Uptake.ratio.gg[,5:7], 1, mean)

Uptake.ratio.gg$sd_long<- apply(Uptake.ratio.gg[,2:4], 1, sd)

Uptake.ratio.gg$sd_short<- apply(Uptake.ratio.gg[,5:7], 1, sd)



#####################    make a flag to indicate low input positions    #########################


#get a  mean of input samples, CAREFULL column numbers might change by function fread

raw.depth.inputs.gg$short.mean<-  apply(raw.depth.inputs.gg[,5:7], 1, mean)  
raw.depth.inputs.gg$large.mean<-  apply(raw.depth.inputs.gg[,2:4], 1, mean)



flag_small<-  raw.depth.inputs.gg$short.mean[1:length(raw.depth.inputs.gg$short.mean)]>=10  #make the flag vector of low input positions
flag_large<-  raw.depth.inputs.gg$large.mean[1:length(raw.depth.inputs.gg$large.mean)]>=10  #make the flag vector of low input positions

Uptake.ratio.gg$flag_small <- flag_small
Uptake.ratio.gg$flag_large <- flag_large

################################      make a normalized depth dataframe and a pseudocount uptake ratio    #########################################
# make a normalized depth dataframe
norm.ratios.gg<-  data.table(pos = remix14.UP04$pos, UP04 = remix14.UP04$UP04.n, UP05 = remix14.UP05$UP05.n, UP06 = remix14.UP06$UP06.n, UP10 = remix16.UP10$UP10.n, UP11 = remix16.UP11$UP11.n, UP12 = remix16.UP12$UP12.n, UP14.UP04 = remix14.UP04$input.n, UP14.UP05 = remix14.UP05$input.n, UP14.UP06 = remix14.UP06$input.n, UP16.UP10 = remix16.UP10$input.n, UP16.UP11 = remix16.UP11$input.n, UP16.UP12 = remix16.UP12$input.n) 
# save
write.csv(norm.ratios.gg, "~/DNA_uptake/datasets/tables/norm.ratios.gg.corrected.csv")
norm.ratios.gg<- read.csv("~/DNA_uptake/datasets/tables/norm.ratios.gg.csv")

m<- apply(norm.ratios.gg[,2:13], MARGIN = 2, mean)


# Add pseudocounts to normalized depth will help removing inout 0 positions. Input 0 positions will affect the peak finder  analysis. 
str(norm.ratios.gg)
which(norm.ratios.gg$UP16.UP10 == 0) #check if there are 0s in the input
pos<- norm.ratios.gg$pos 
UP10.ratio<- norm.ratios.gg$UP10/(norm.ratios.gg$UP16.UP10 + 1)
UP11.ratio<- norm.ratios.gg$UP11/(norm.ratios.gg$UP16.UP11 + 1)
UP12.ratio<- norm.ratios.gg$UP12/(norm.ratios.gg$UP16.UP12 + 1)

pseudo.ratios.gg<- data.frame( pos = pos, UP10.ratio = UP10.ratio, UP11.ratio = UP11.ratio, UP12.ratio = UP12.ratio  ) #make data.frame

pseudo.ratios.gg$ratio_short<- apply(pseudo.ratios.gg[,2:4], FUN = mean, MARGIN = 1)

str(pseudo.ratios.gg)

write.csv(pseudo.ratios.gg, "~/DNA_uptake/datasets/tables/pseudo.ratios.gg.csv")


#############################save or load dataframes (if needed)##############################

str(raw.depth.samples.gg)
str(raw.depth.inputs.gg)
str(Uptake.ratio.gg)


#save dataframes
write.csv(raw.depth.samples.gg, "~/DNA_uptake/datasets/tables/raw.depth.samples.gg.csv")
write.csv(raw.depth.inputs.gg, "~/DNA_uptake/datasets/tables/raw.depth.inputs.gg.csv")
write.csv(Uptake.ratio.gg, "~/DNA_uptake/datasets/tables/Uptake.ratio.gg.corrected.csv")

# load dataframes (OPTIONAL CODE)
raw.depth.samples.gg<- fread("~/DNA_uptake/datasets/tables/raw.depth.samples.gg.csv") #load new input samples fro UP01   
raw.depth.inputs.gg<- fread("~/DNA_uptake/datasets/tables/raw.depth.inputs.gg.csv") #load new input samples fro UP01   
Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv") #load new input samples fro UP01   

dif<- Uptake.ratio.gg.c$ratio_short - Uptake.ratio.gg$ratio_short

summary(dif)


# remove the extra column introduced by fread
Uptake.ratio.gg$V1 <- NULL
raw.depth.inputs.gg$V1 <- NULL
raw.depth.samples.gg$V1 <- NULL

##############################################################################################

#####################################      Some Pretty uptake maps for short fragments   ###########################

# save path where figures will be saved
### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/Uptake_maps/"
# Set directory
setwd(whereami)



u<- Uptake.ratio.gg %>%  dplyr::slice(10000:20000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_line(aes(colour = flag_small))+
  geom_ribbon(aes(ymin=ratio_short-sd_short, ymax=ratio_short+sd_short), color = "blue", alpha =0.3) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(10000 , 20000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("short.ratio.gg.sd2","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.gg %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.gg","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  


summary(Uptake.ratio.gg$ratio_short)


  u<- Uptake.ratio.gg %>%  
    ggplot(aes(x = pos,y = ratio_short)) +
    geom_point(shape = 20, size = 1, aes(colour = flag_small))+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
    scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
    scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16,24,48, 96), limits = c(NA, 96), expand = c(0, 0))+
    labs(x = "PittGG genome positions", y = "uptake ratio") +
    ggtitle(" Uptake ratio short donor PittGG DNA") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "bottom",
          panel.grid.minor = element_line(colour="white", size=0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text  = element_text(size=18),
          axis.title = element_text(size = 18, face = "bold")) 
  
  file_name = paste("short.ratio.gg","whole_genome", "log", "2", "tiff", sep=".")
  tiff(file_name, width = 1000, height = 700, units = "px")
  print(u)
  dev.off()  
    
u<- Uptake.ratio.gg %>%  
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0,80 ), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.gg","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(u)
dev.off()  



u<- Uptake.ratio.gg %>%  dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.gg","100kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()



u<- Uptake.ratio.gg %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16), limits = c(NA, 12), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio short donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.gg.log","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  


u<- Uptake.ratio.gg %>%  dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16,24,48), limits = c(NA, 48), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio short donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.gg.log.3","100kb","2", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  





u<- Uptake.ratio.gg  %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("short.ratio.gg","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(u)
dev.off()  



#####################################      Some Pretty uptake maps for large fragments   ###########################

# save path where figures will be saved
### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/Uptake_maps/"
# Set directory
setwd(whereami)



u<- Uptake.ratio.gg %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_line(aes(colour = flag_large))+
  geom_errorbar(aes(ymin=ratio_long-sd_long, ymax=ratio_long+sd_long), colour = "grey", alpha = 1) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("large.ratio.gg.errorbars","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.gg %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("large.ratio.gg","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.gg %>%  dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("large.ratio.gg","100kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()



u<- Uptake.ratio.gg %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16), limits = c(NA, 12), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio large donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("large.ratio.gg.log","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  




u<- Uptake.ratio.gg  %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("large.ratio.gg","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.gg  %>%  
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
  scale_y_continuous(trans = "log2", limits = c(0.0001220703, 64), breaks = c(0.0001220703, 0.001953125,0.0625,2,64), labels = c(0.0001, 0.0020,0.0625,2,64), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio small donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("./small.ratio.gg.log2","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1200, height = 700, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.gg  %>%  
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1887050, 200000), expand = c(0, 0))+
  scale_y_continuous(trans = "log2", limits = c(0.0001220703, 64), breaks = c(0.0001220703, 0.001953125,0.0625,2,64), labels = c(0.0001, 0.0020,0.0625,2,64), expand = c(0, 0))+
  labs(x = "PittGG genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio long donor PittGG DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("./long.ratio.gg.log2","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1200, height = 700, units = "px")
print(u)
dev.off()


#######################################################################

#############      compare new vs old ratios    #######################


Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv") #load new input samples fro UP01   

Uptake.ratio.gg.c<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.corrected.csv") #load new input samples fro UP01   


fold<- Uptake.ratio.gg.c$ratio_short/Uptake.ratio.gg$ratio_short

test<- data.frame(method1_ratio = Uptake.ratio.gg$ratio_short, fold_change = fold)


summary(test$method1_ratio)


u<- ggplot() +
  geom_point(aes(x = method1_ratio,y = fold_change), shape = ".", size = 1, colour = "blue", data = test)+
  scale_x_continuous(limits = c(0, 40),breaks = seq(0 , 40, 5), expand = c(0, 0))+
  labs(x = "PittGG method 1 uptake ratio", y = " fold change method 2/method 1") +
  ggtitle(" PittGG fold difference method 2 vs method 1") +
  theme_grey()+
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("./method.diff","whole_genome","old vs new", "tiff", sep=".")
tiff(file_name, width = 1200, height = 700, units = "px")
print(u)
dev.off()


##########################################   add row input normalized cov  #######################################

unique(samples.gg$sample)


up14<- dplyr::filter(samples.gg, sample == "UP14") #subset UP01 

str(up14)

up14.mapped<- sum(up14$depth)/length(up14$depth)  
up14.n<- up14$depth/up14.mapped



up16<- dplyr::filter(samples.gg, sample == "UP16") #subset UP01 

str(up16)

up16.mapped<- sum(up16$depth)/length(up16$depth)  
up16.n<- up16$depth/up16.mapped

mean(up16.n)


str(norm.ratios.gg)


norm.ratios.gg$UP14<- up14.n

norm.ratios.gg$UP16 <- up16.n

write.csv(norm.ratios.gg, "./datasets/new_norm_datasets/norm.ratios.gg.corrected.csv")

