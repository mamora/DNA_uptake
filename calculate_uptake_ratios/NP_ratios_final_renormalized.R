

#####################################    Goal     ######################################### 

#To generate uptake ratios for small and large fragment data for NP donor DNA samples




#############################load packages##############################
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
library(tidyr)
#########################################################################

#############################load dataframes##############################
remix13.UP01<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX13.UP01.I.np.depth.bed") #load new input samples fro UP01   
colnames(remix13.UP01)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
remix13.UP01$dif<- remix13.UP01$depth.1/((remix13.UP01$depth.2+remix13.UP01$depth.3+remix13.UP01$depth.4)/3)
samples.np<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/samples.all.np.bed.csv") #load all samples from NP donor DNA, input included
colnames(samples.np)<- c("sample","genome","pos","depth") # add headers
up01<- dplyr::filter(samples.np, sample == "UP01") #subset UP01 
remix13.UP01$UP01<- up01$depth
remix13.UP01$input.mean<- (remix13.UP01$depth.1+remix13.UP01$depth.2+remix13.UP01$depth.3+remix13.UP01$depth.4)/4
#########################################################################

n<- which(remix13.UP01$input.mean <= 10)

(74/length(remix13.UP01$input.mean)*100)

#############################check dataframes##############################
str(samples.np)
str(remix13.UP01)
unique(samples.np$sample) #check sample names
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization#### ##########################
up13.up01.mapped<- sum(remix13.UP01$input.mean)/ length(remix13.UP01$input.mean)
up01.all.mapped<- sum(remix13.UP01$UP01)/length(remix13.UP01$input.mean)
remix13.UP01$UP01.n<- (remix13.UP01$UP01)/up01.all.mapped
remix13.UP01$input.n<- (remix13.UP01$input.mean)/up13.up01.mapped
remix13.UP01$UP01.ratio<- remix13.UP01$UP01.n/remix13.UP01$input.n
mean(remix13.UP01$UP01.ratio)
mean(remix13.UP01$UP01.n)
#########################################################################

mean(remix13.UP01$UP01.n)
mean(remix13.UP01$input.n)


##############################plot coverage rdnp all first 1kb all samples KW20 #####################################

theme_set(theme_grey()) #set theme

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/samples_figs/"
# Set directory
setwd(whereami)


u<- remix13.UP01 %>%  dplyr::slice(1:10000) %>% 
    ggplot(aes(x = pos,y = UP01.ratio)) +
  geom_point(shape = 20, size = 1)+
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP01") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("UP01.ratio","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  


u<- remix13.UP01 %>% 
  ggplot(aes(x = pos,y = UP01.ratio)) +
  geom_point(shape = 20, size = 1)+
  scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP01") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))
file_name = paste("UP01.ratio","whole_genome", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off() 


###################################     UP02       ######################################

#############################load dataframes##############################
remix13.UP02<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX13.UP02.I.np.depth.bed") #load new input samples for UP02   
colnames(remix13.UP02)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up02<- dplyr::filter(samples.np, sample == "UP02") #subset UP01 
remix13.UP02$UP02<- up02$depth
remix13.UP02$input.mean<- (remix13.UP02$depth.1+remix13.UP02$depth.2+remix13.UP02$depth.3+remix13.UP02$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(samples.np)
str(remix13.UP02)
unique(samples.np$sample) #check sample names
#########################################################################

sum(remix13.UP01$input.mean)

##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up13.up02.mapped<- sum(remix13.UP02$input.mean)/length(remix13.UP01$input.mean) # total mapped reads, from sambamba flagstats files 
up02.all.mapped<- sum(remix13.UP02$UP02)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix13.UP02$UP02.n<- (remix13.UP02$UP02)/up02.all.mapped
remix13.UP02$input.n<- (remix13.UP02$input.mean)/up13.up02.mapped
remix13.UP02$UP02.ratio<- remix13.UP02$UP02.n/remix13.UP02$input.n
#########################################################################


sum(remix13.UP02$UP02.n)

mean(remix13.UP02$input.n)


#################################################################################################

###################################     UP03       ######################################

#############################load dataframes##############################
remix13.UP03<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX13.UP03.I.np.depth.bed") #load new input samples for UP03   
colnames(remix13.UP03)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up03<- dplyr::filter(samples.np, sample == "UP03") #subset UP03 
remix13.UP03$UP03<- up03$depth
remix13.UP03$input.mean<- (remix13.UP03$depth.1+remix13.UP03$depth.2+remix13.UP03$depth.3+remix13.UP03$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(samples.np)
str(remix13.UP03)
unique(samples.np$sample) #check sample names
#########################################################################



##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up13.up03.mapped<- sum(remix13.UP03$input.mean)/length(remix13.UP01$input.mean) # total mapped reads, from sambamba flagstats files 
up03.all.mapped<- sum(remix13.UP03$UP03)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix13.UP03$UP03.n<- (remix13.UP03$UP03)/up03.all.mapped
remix13.UP03$input.n<- (remix13.UP03$input.mean)/up13.up03.mapped
remix13.UP03$UP03.ratio<- remix13.UP03$UP03.n/remix13.UP03$input.n
#########################################################################

####################################     save in dataframes     ###############################
raw.depth.samples<- data.table(pos = remix13.UP01$pos, UP01 = remix13.UP01$UP01, UP02 = remix13.UP02$UP02, UP03 = remix13.UP03$UP03) 
raw.depth.inputs<-  data.table(pos = remix13.UP01$pos, UP13.UP01 = remix13.UP01$input.mean, UP13.UP02 = remix13.UP02$input.mean, UP13.UP03 = remix13.UP03$input.mean) 
Uptake.ratio.np<-  data.table(pos = remix13.UP01$pos, UP01.ratio = remix13.UP01$UP01.ratio, UP02.ratio = remix13.UP02$UP02.ratio, UP03.ratio = remix13.UP03$UP03.ratio) 




#################################      UP07       ############################################

#############################load dataframes##############################
remix15.UP07<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX15.UP07.I.np.depth.bed") #load new input samples fro UP07   
colnames(remix15.UP07)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up07<- dplyr::filter(samples.np, sample =="UP07") #subset UP07 
up15<- dplyr::filter(samples.np, sample =="UP15") #subset UP07 
remix15.UP07$dif<- remix15.UP07$depth.1/((remix15.UP07$depth.2+remix15.UP07$depth.3+remix15.UP07$depth.4)/3)
remix15.UP07$UP07<- up07$depth
remix15.UP07$input.mean<- (remix15.UP07$depth.1+remix15.UP07$depth.2+remix15.UP07$depth.3+remix15.UP07$depth.4)/4

#########################################################################


#############################check dataframes##############################
str(remix15.UP07)
str(samples.np)
str(up15)
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up15.up07.mapped<- sum(remix15.UP07$input.mean)/length(remix13.UP01$input.mean) # total mapped reads, from sambamba flagstats files 
up07.all.mapped<- sum(remix15.UP07$UP07)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix15.UP07$UP07.n<- (remix15.UP07$UP07)/up07.all.mapped
remix15.UP07$input.n<- (remix15.UP07$input.mean)/up15.up07.mapped
remix15.UP07$UP07.ratio<- remix15.UP07$UP07.n/remix15.UP07$input.n
#########################################################################

sum(remix15.UP07$UP07)
sum(remix15.UP07$input.mean)

# save in dataframes
raw.depth.samples$UP07<- remix15.UP07$UP07
raw.depth.inputs$UP15.UP07<- remix15.UP07$input.mean
Uptake.ratio.np$UP07.ratio<- remix15.UP07$UP07.ratio


##############################plot coverage rdnp all first 1kb all samples KW20 #####################################

theme_set(theme_grey()) #set theme

### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/samples_figs/"
# Set directory
setwd(whereami)


u<- remix15.UP07 %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = UP07.ratio)) +
  geom_point(shape = 20, size = 1)+
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP07") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("UP07.ratio.2","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, bg = "transparent", units = "px")
print(u)
dev.off() 



u<- remix15.UP07 %>%  dplyr::slice(1:10000) %>%
  ggplot(aes(x = pos,y = UP07.ratio)) +
  geom_point(shape = 20, size = 1)+
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16), limits = c(NA, 12))+    
  labs(x = "NP genome positions", y = "log uptake ratio") +
  ggtitle(" Uptake ratio UP07 log") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("UP07.ratio.log","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, bg = "transparent", units = "px")
print(u)
dev.off()  


u<- remix15.UP07 %>% 
  ggplot(aes(x = pos,y = UP07.ratio)) +
  geom_point(shape = ".", size = 1)+
  scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio UP07") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))
file_name = paste("UP07.ratio","whole_genome", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off() 

u<- uptake  %>% dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio)) +
  geom_point(shape = 20, size = 1)+
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio old data") +
  theme()+
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
      panel.border = element_blank(),
      legend.key = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA))
file_name = paste("old.ratio.np.t","10kb", "tiff", sep=".")
png(file_name, width = 800, height = 500,bg = "transparent", units = "px")
print(u)
dev.off()  

  
  u<- uptake  %>% dplyr::slice(1:10000)
  u$ratio.2<- u$ratio*2
  ggplot(aes(x = pos,y = ratio.2), data = u) +
  geom_point(shape = 20, size = 1, color = "blue")+
  geom_point(aes(x = pos, y = UP07.ratio), color = "red", shape = 20, size = 1, data = u1)+
    scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16), limits = c(NA, 12))+
      labs(x = "NP genome positions", y = "uptake ratio log scale") +
  ggtitle(" Uptake ratio old (blue) and new (red) data log scale") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("old.new.ratio.np.log*2","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off() 
  
############################################################################################

######################################      UP08        ###################################################



#############################load dataframes##############################
remix15.UP08<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX15.UP08.I.np.depth.bed") #load new input samples fro UP08  
colnames(remix15.UP08)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up08<- dplyr::filter(samples.np, sample =="UP08") #subset UP07 
remix15.UP08$UP08<- up08$depth
remix15.UP08$input.mean<- (remix15.UP08$depth.1+remix15.UP08$depth.2+remix15.UP08$depth.3+remix15.UP08$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(remix15.UP08)
str(samples.np)
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up15.up08.mapped<- sum(remix15.UP08$input.mean)/length(remix13.UP01$input.mean)  # total mapped reads, from sambamba flagstats files 
up08.all.mapped<- sum(remix15.UP08$UP08)/length(remix13.UP01$input.mean)# total mapped reads from summary table
remix15.UP08$UP08.n<- (remix15.UP08$UP08)/up08.all.mapped
remix15.UP08$input.n<- (remix15.UP08$input.mean)/up15.up08.mapped
remix15.UP08$UP08.ratio<- remix15.UP08$UP08.n/remix15.UP08$input.n
#########################################################################

sum(remix15.UP08$UP08)
sum(remix15.UP08$input.mean)


raw.depth.samples$UP08<- remix15.UP08$UP08
raw.depth.inputs$UP15.UP08<- remix15.UP08$input.mean
Uptake.ratio.np$UP08.ratio<- remix15.UP08$UP08.ratio

#######################################################################################################


######################################      UP09        ###################################################



#############################load dataframes##############################
remix15.UP09<- fread("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/REMIX/REMIX_files/all/REMIX15.UP09.I.np.depth.bed") #load new input samples fro UP09   
colnames(remix15.UP09)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up09<- dplyr::filter(samples.np, sample =="UP09") #subset UP07 
remix15.UP09$UP09<- up09$depth
remix15.UP09$input.mean<- (remix15.UP09$depth.1+remix15.UP09$depth.2+remix15.UP09$depth.3+remix15.UP09$depth.4)/4
#########################################################################


#############################check dataframes##############################
str(remix15.UP09)
str(samples.np)
#########################################################################


##Normalize the raw reads to depth per millon reads##
#############################normalization##############################
up15.up09.mapped<- sum(remix15.UP09$input.mean)/length(remix13.UP01$input.mean)  # total mapped reads, from sambamba flagstats files 
up09.all.mapped<- sum(remix15.UP09$UP09)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix15.UP09$UP09.n<- (remix15.UP09$UP09)/up09.all.mapped
remix15.UP09$input.n<- (remix15.UP09$input.mean)/up15.up09.mapped
remix15.UP09$UP09.ratio<- remix15.UP09$UP09.n/remix15.UP09$input.n
#########################################################################

sum(remix15.UP09$UP09)
sum(remix15.UP09$input.mean)


raw.depth.samples$UP09<- remix15.UP09$UP09
raw.depth.inputs$UP15.UP09<- remix15.UP09$input.mean
Uptake.ratio.np$UP09.ratio<- remix15.UP09$UP09.ratio

#######################################################################################################



##################calculate mean and standard deviation of the three replicates#######################

Uptake.ratio.np$ratio_long<- apply(Uptake.ratio.np[,2:4], 1, mean)

Uptake.ratio.np$ratio_short<- apply(Uptake.ratio.np[,5:7], 1, mean)

Uptake.ratio.np$sd_long<- apply(Uptake.ratio.np[,2:4], 1, sd)

Uptake.ratio.np$sd_short<- apply(Uptake.ratio.np[,5:7], 1, sd)

m1<- apply(Uptake.ratio.np[,2:9], MARGIN = 2, mean)

#####################    make a flag to indicate low input positions    #########################


#get a  mean of input samples, CAREFULL column numbers might change by function fread

raw.depth.inputs$short.mean<-  apply(raw.depth.inputs[,5:7], 1, mean)  
raw.depth.inputs$large.mean<-  apply(raw.depth.inputs[,2:4], 1, mean)



flag_small<-  raw.depth.inputs$short.mean[1:length(raw.depth.inputs$short.mean)]>=10  #make the flag vector of low input positions
flag_large<-  raw.depth.inputs$large.mean[1:length(raw.depth.inputs$large.mean)]>=10  #make the flag vector of low input positions

Uptake.ratio.np$flag_small <- flag_small
Uptake.ratio.np$flag_large <- flag_large

################################      make a normalized depth dataframe and a pseudocount uptake ratio    #########################################
# make a normalized depth dataframe
norm.ratios.np<-  data.table(pos = remix13.UP01$pos, UP01 = remix13.UP01$UP01.n, UP02 = remix13.UP02$UP02.n, UP03 = remix13.UP03$UP03.n, UP07 = remix15.UP07$UP07.n, UP08 = remix15.UP08$UP08.n, UP09 = remix15.UP09$UP09.n, UP13.UP01 = remix13.UP01$input.n, UP13.UP02 = remix13.UP02$input.n, UP13.UP03 = remix13.UP03$input.n, UP15.UP07 = remix15.UP07$input.n, UP15.UP08 = remix15.UP08$input.n, UP15.UP09 = remix15.UP09$input.n) 
# save
write.csv(norm.ratios.np, "~/DNA_uptake/datasets/tables/norm.ratios.np.corrected.csv")
# mean
m<- apply(norm.ratios.np[,2:12], MARGIN = 2, mean)


# Add pseudocounts to normalized depth will help removing inout 0 positions. Input 0 positions will affect the peak finder  analysis. 
str(norm.ratios.np)
pseudo.ratios.np<- data.frame() #make empty data.frame
which(norm.ratios.np$UP15.UP09 < 1) #check if there are 0s in the input
norm.ratios.np$UP15.UP09[681619]
pseudo.ratios.np$UP01.ratio<- norm.ratios.np$UP01/(norm.ratios.np$UP13.UP01 + 1)


#THIS SECTION IS SUSPENDED SINCE THERE SEEMS TO BE NO INPUT 0 POSIIONS, SO I WON'T NEED TO ADD A PSEUDOCOUNT 

#############################save or load dataframes (if needed)##############################

str(raw.depth.samples)
str(raw.depth.inputs)
str(Uptake.ratio.np)


#save dataframes
write.csv(raw.depth.samples, "~/DNA_uptake/datasets/tables/raw.depth.samples.csv")
write.csv(raw.depth.inputs, "~/DNA_uptake/datasets/tables/raw.depth.inputs.csv")
write.csv(Uptake.ratio.np, "~/DNA_uptake/datasets/tables/Uptake.ratio.np.corrected.csv")

# load dataframes (OPTIONAL CODE)
raw.depth.samples<- fread("~/DNA_uptake/datasets/tables/raw.depth.samples.csv") #load new uptake samples   
raw.depth.inputs<- fread("~/DNA_uptake/datasets/tables/raw.depth.inputs.csv") #load new input samples  
Uptake.ratio.np<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.np.csv") #load uptake ratios   
Uptake.ratio.gg<- fread("~/DNA_uptake/datasets/tables/Uptake.ratio.gg.csv")    



# remove the extra column introduced by fread
Uptake.ratio.np$V1 <- NULL
raw.depth.inputs$V1 <- NULL
raw.depth.samples$V1 <- NULL

# function for normalizing predicted uptake
norm<- function (data = data){
  h_pe<- mean(data, na.rm = FALSE)
  s_pe<- (data * 1)/h_pe
  return(s_pe)
}

# normalize predicted uptake to a mean of 1
Uptake.ratio.np$ratio_short<- norm(data = Uptake.ratio.np$ratio_short)
Uptake.ratio.np$ratio_long<- norm(data = Uptake.ratio.np$ratio_long)


mean(Uptake.ratio.np$ratio_short)
mean(Uptake.ratio.np$ratio_long)

#read uptake file
write.csv(Uptake.ratio.np, "./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")      

##############################################################################################

#####################################      Some Pretty uptake maps   ###########################

# save path where figures will be saved
### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/Uptake_maps/"
# Set directory
setwd(whereami)



u<- Uptake.ratio.np %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_line()+
  geom_ribbon(aes(ymin=ratio_short-sd_short, ymax=ratio_short+sd_short), color = "blue", alpha =0.3) +
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor NP DNA with standard deviation") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("short.ratio.np.sd","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.np %>%  dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_line()+
  geom_ribbon(aes(ymin=ratio_short-sd_short, ymax=ratio_short+sd_short), color = "blue", alpha =0.3) +
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor NP DNA with standard deviation") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("short.ratio.np.sd","100kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()

u<- Uptake.ratio.np %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.np","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.np %>%  dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.np","100kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()

#join gg and np dataframes to force ggplot to use the same y axis
joined<- data.frame(pos = Uptake.ratio.np$pos[1:100000], ratio.np = Uptake.ratio.np$ratio_short[1:100000], ratio.gg = Uptake.ratio.gg$ratio_short[1:100000])  

joined.tidy<- joined %>% tidyr::gather(sample, ratio, -pos)



u<- ggplot(aes(x = pos,y = ratio), data = joined.tidy) +
  geom_point(shape = 20, size = 1, colour = "blue")+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16, 24,48), limits = c(NA, 48), expand = c(0.0005, 0))+
  facet_grid(sample ~.) +
  labs(x = "genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.np.gg.log","100kb", "tiff", sep=".")
tiff(file_name, width = 1400, height = 1000, units = "px")
print(u)
dev.off() 



  
  u<- Uptake.ratio.np %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16, 24,48), limits = c(NA, 12), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.np.log","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  
  

u<- Uptake.ratio.np %>%  dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(trans = log2_trans(), breaks = c(NA,0.0005,0.005,0.05,0.5,1,2,4,8,16,24,48), limits = c(NA, 48), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("short.ratio.np.log","100kb","2", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()
  
  
  
u<- Uptake.ratio.np  %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio short donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("short.ratio.np","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(u)
dev.off()  


#####################################      Some Pretty uptake maps for large fragments   ###########################

# save path where figures will be saved
### YOUR PATH IS DIFFERENT THAN THIS ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/Uptake_maps/"
# Set directory
setwd(whereami)



u<- Uptake.ratio.np %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_line()+
  geom_ribbon(aes(ymin=ratio_long-sd_long, ymax=ratio_long+sd_long), color = "blue", alpha =0.3) +
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("large.ratio.np.sd","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.np %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("large.ratio.np","10kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.np %>%  dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = 20, size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("large.ratio.np","100kb", "tiff", sep=".")
tiff(file_name, width = 800, height = 500, units = "px")
print(u)
dev.off()


log2(128)

0.0078125/2

0.00390625/2

0.0002441406/2


u<- Uptake.ratio.np  %>% dplyr::slice(1:100000) %>% 
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 100000, 10000), expand = c(0, 0))+
  scale_y_continuous(trans = "log2", limits = c(0.0001220703, 64), breaks = c(0.0001220703, 0.001953125,0.0625,2,64), labels = c(0.0001, 0.0020,0.0625,2,64), expand = c(0, 0))+
    labs(x = "NP genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio large donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("./small.ratio.np.log","test","log2", "tiff", sep=".")
tiff(file_name, width = 800, height = 800, units = "px")
print(u)
dev.off()  


u<- Uptake.ratio.np  %>%  
  ggplot(aes(x = pos,y = ratio_short)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_small))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
  scale_y_continuous(trans = "log2", limits = c(0.0001220703, 64), breaks = c(0.0001220703, 0.001953125,0.0625,2,64), labels = c(0.0001, 0.0020,0.0625,2,64), expand = c(0, 0))+
  labs(x = "NP genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio small donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("./small.ratio.np.log2","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1200, height = 700, units = "px")
print(u)
dev.off()  

u<- Uptake.ratio.np  %>%  
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
  scale_y_continuous(trans = "log2", limits = c(0.0001220703, 64), breaks = c(0.0001220703, 0.001953125,0.0625,2,64), labels = c(0.0001, 0.0020,0.0625,2,64), expand = c(0, 0))+
  labs(x = "NP genome positions", y = " log2 of uptake ratio") +
  ggtitle(" Uptake ratio long donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 

file_name = paste("./long.ratio.np.log2","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1200, height = 700, units = "px")
print(u)
dev.off()

u<- Uptake.ratio.np  %>% 
  ggplot(aes(x = pos,y = ratio_long)) +
  geom_point(shape = ".", size = 1, aes(colour = flag_large))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Input coverage", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("TRUE" = "more than 10" , "FALSE" = "less or equal to 10")) + 
  scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake ratio") +
  ggtitle(" Uptake ratio large donor NP DNA") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("large.ratio.np","whole_genome", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(u)
dev.off()  






