################################################################################################
#     This script takes the raw bed files from the sequence handling scripts and calculates    #
#     normalized uptake ratios and well as normalized depth files for NP genome                #
################################################################################################

# load packages
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
library(tidyr)

# The way this script will work is that it will take a sample bed file (see table S1) and it will:
# 1. calculate the depth of donor sample, 2. calculate the depth of input sample, 3. normalize both
# 4. divide donor normalize depth by input normalize depth to generate an uptake ratio.

# lastly, uptake ratios of all mean replicates will be generated and will be normalized 
# to a mean of 1

##################################################################################################
#   The following lines will calculate uptake ratio for each sample at a time                    #
##################################################################################################


###################################     UP01       ######################################

# load input and donor bed files and calculate uptake for UP01
remix13.UP01<- fread("./REMIX/REMIX_files/all/REMIX13.UP01.I.np.depth.bed") #load new input samples fro UP01   
colnames(remix13.UP01)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
remix13.UP01$dif<- remix13.UP01$depth.1/((remix13.UP01$depth.2+remix13.UP01$depth.3+remix13.UP01$depth.4)/3)
samples.np<- fread("./REMIX/REMIX_files/all/samples.all.np.bed.csv") #load all samples from NP donor DNA, input included
colnames(samples.np)<- c("sample","genome","pos","depth") # add headers
up01<- dplyr::filter(samples.np, sample == "UP01") #subset UP01 
remix13.UP01$UP01<- up01$depth
remix13.UP01$input.mean<- (remix13.UP01$depth.1+remix13.UP01$depth.2+remix13.UP01$depth.3+remix13.UP01$depth.4)/4

# check dataframes
str(samples.np)
str(remix13.UP01)
unique(samples.np$sample) #check sample names

# Normalize the raw reads to depth per millon reads##
up13.up01.mapped<- sum(remix13.UP01$input.mean)/ length(remix13.UP01$input.mean)
up01.all.mapped<- sum(remix13.UP01$UP01)/length(remix13.UP01$input.mean)
remix13.UP01$UP01.n<- (remix13.UP01$UP01)/up01.all.mapped
remix13.UP01$input.n<- (remix13.UP01$input.mean)/up13.up01.mapped
remix13.UP01$UP01.ratio<- remix13.UP01$UP01.n/remix13.UP01$input.n



###################################     UP02       ######################################

# load input and donor bed files and calculate uptake for UP02
remix13.UP02<- fread("./REMIX/REMIX_files/all/REMIX13.UP02.I.np.depth.bed") #load new input for UP02  
# add headers
colnames(remix13.UP02)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") 
up02<- dplyr::filter(samples.np, sample == "UP02") #subset donor UP02 
remix13.UP02$UP02<- up02$depth
remix13.UP02$input.mean<- (remix13.UP02$depth.1+remix13.UP02$depth.2+remix13.UP02$depth.3+remix13.UP02$depth.4)/4

# check dataframes
str(samples.np)
str(remix13.UP02)
unique(samples.np$sample) #check sample names

# Normalize the raw reads to depth per millon reads##
up13.up02.mapped<- sum(remix13.UP02$input.mean)/length(remix13.UP01$input.mean) # total mapped reads, from sambamba flagstats files 
up02.all.mapped<- sum(remix13.UP02$UP02)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix13.UP02$UP02.n<- (remix13.UP02$UP02)/up02.all.mapped
remix13.UP02$input.n<- (remix13.UP02$input.mean)/up13.up02.mapped
remix13.UP02$UP02.ratio<- remix13.UP02$UP02.n/remix13.UP02$input.n


###################################     UP03       ######################################

# load input and donor bed files and calculate uptake for UP03
remix13.UP03<- fread("./REMIX/REMIX_files/all/REMIX13.UP03.I.np.depth.bed") #load new input samples for UP03   
# add headers
colnames(remix13.UP03)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") 
up03<- dplyr::filter(samples.np, sample == "UP03") #subset UP03 
remix13.UP03$UP03<- up03$depth
remix13.UP03$input.mean<- (remix13.UP03$depth.1+remix13.UP03$depth.2+remix13.UP03$depth.3+remix13.UP03$depth.4)/4


# check dataframes
str(samples.np)
str(remix13.UP03)
unique(samples.np$sample) #check sample names

# Normalize the raw reads to depth per millon reads##
up13.up03.mapped<- sum(remix13.UP03$input.mean)/length(remix13.UP01$input.mean) # total mapped reads, from sambamba flagstats files 
up03.all.mapped<- sum(remix13.UP03$UP03)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix13.UP03$UP03.n<- (remix13.UP03$UP03)/up03.all.mapped
remix13.UP03$input.n<- (remix13.UP03$input.mean)/up13.up03.mapped
remix13.UP03$UP03.ratio<- remix13.UP03$UP03.n/remix13.UP03$input.n

####################################     save in dataframes     ###############################
raw.depth.samples<- data.table(pos = remix13.UP01$pos, UP01 = remix13.UP01$UP01, UP02 = remix13.UP02$UP02, UP03 = remix13.UP03$UP03) 
raw.depth.inputs<-  data.table(pos = remix13.UP01$pos, UP13.UP01 = remix13.UP01$input.mean, UP13.UP02 = remix13.UP02$input.mean, UP13.UP03 = remix13.UP03$input.mean) 
Uptake.ratio.np<-  data.table(pos = remix13.UP01$pos, UP01.ratio = remix13.UP01$UP01.ratio, UP02.ratio = remix13.UP02$UP02.ratio, UP03.ratio = remix13.UP03$UP03.ratio) 

#################################      UP07       ############################################

# load input and donor bed files and calculate uptake for UP07
remix15.UP07<- fread("./REMIX/REMIX_files/all/REMIX15.UP07.I.np.depth.bed") #load new input samples for UP07   
colnames(remix15.UP07)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up07<- dplyr::filter(samples.np, sample =="UP07") #subset UP07 
remix15.UP07$dif<- remix15.UP07$depth.1/((remix15.UP07$depth.2+remix15.UP07$depth.3+remix15.UP07$depth.4)/3)
remix15.UP07$UP07<- up07$depth
remix15.UP07$input.mean<- (remix15.UP07$depth.1+remix15.UP07$depth.2+remix15.UP07$depth.3+remix15.UP07$depth.4)/4

# check dataframes
str(remix15.UP07)
str(samples.np)
str(up15)

# Normalize the raw reads to depth per millon reads##
up15.up07.mapped<- sum(remix15.UP07$input.mean)/length(remix13.UP01$input.mean) # total mapped reads, from sambamba flagstats files 
up07.all.mapped<- sum(remix15.UP07$UP07)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix15.UP07$UP07.n<- (remix15.UP07$UP07)/up07.all.mapped
remix15.UP07$input.n<- (remix15.UP07$input.mean)/up15.up07.mapped
remix15.UP07$UP07.ratio<- remix15.UP07$UP07.n/remix15.UP07$input.n

# save in dataframes
raw.depth.samples$UP07<- remix15.UP07$UP07
raw.depth.inputs$UP15.UP07<- remix15.UP07$input.mean
Uptake.ratio.np$UP07.ratio<- remix15.UP07$UP07.ratio

######################################      UP08        ###################################################

# load input and donor bed files and calculate uptake for UP08
remix15.UP08<- fread("./REMIX/REMIX_files/all/REMIX15.UP08.I.np.depth.bed") #load new input samples fro UP08  
colnames(remix15.UP08)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up08<- dplyr::filter(samples.np, sample =="UP08") #subset UP07 
remix15.UP08$UP08<- up08$depth
remix15.UP08$input.mean<- (remix15.UP08$depth.1+remix15.UP08$depth.2+remix15.UP08$depth.3+remix15.UP08$depth.4)/4
#########################################################################


# check dataframes
str(remix15.UP08)
str(samples.np)


# Normalize the raw reads to depth per millon reads##
up15.up08.mapped<- sum(remix15.UP08$input.mean)/length(remix13.UP01$input.mean)  # total mapped reads, from sambamba flagstats files 
up08.all.mapped<- sum(remix15.UP08$UP08)/length(remix13.UP01$input.mean)# total mapped reads from summary table
remix15.UP08$UP08.n<- (remix15.UP08$UP08)/up08.all.mapped
remix15.UP08$input.n<- (remix15.UP08$input.mean)/up15.up08.mapped
remix15.UP08$UP08.ratio<- remix15.UP08$UP08.n/remix15.UP08$input.n

raw.depth.samples$UP08<- remix15.UP08$UP08
raw.depth.inputs$UP15.UP08<- remix15.UP08$input.mean
Uptake.ratio.np$UP08.ratio<- remix15.UP08$UP08.ratio

######################################      UP09        ###################################################

# load input and donor bed files and calculate uptake for UP09
remix15.UP09<- fread("./REMIX/REMIX_files/all/REMIX15.UP09.I.np.depth.bed") #load new input samples fro UP09   
colnames(remix15.UP09)<- c("genome","start","end","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up09<- dplyr::filter(samples.np, sample =="UP09") #subset UP07 
remix15.UP09$UP09<- up09$depth
remix15.UP09$input.mean<- (remix15.UP09$depth.1+remix15.UP09$depth.2+remix15.UP09$depth.3+remix15.UP09$depth.4)/4

# check dataframes
str(remix15.UP09)
str(samples.np)

# Normalize the raw reads to depth per millon reads##
up15.up09.mapped<- sum(remix15.UP09$input.mean)/length(remix13.UP01$input.mean)  # total mapped reads, from sambamba flagstats files 
up09.all.mapped<- sum(remix15.UP09$UP09)/length(remix13.UP01$input.mean) # total mapped reads from summary table
remix15.UP09$UP09.n<- (remix15.UP09$UP09)/up09.all.mapped
remix15.UP09$input.n<- (remix15.UP09$input.mean)/up15.up09.mapped
remix15.UP09$UP09.ratio<- remix15.UP09$UP09.n/remix15.UP09$input.n

raw.depth.samples$UP09<- remix15.UP09$UP09
raw.depth.inputs$UP15.UP09<- remix15.UP09$input.mean
Uptake.ratio.np$UP09.ratio<- remix15.UP09$UP09.ratio

#save dataframes
write.csv(raw.depth.samples, "./datasets/final_datasets/raw.depth.samples.csv")
write.csv(raw.depth.inputs, "./datasets/final_datasets/raw.depth.inputs.csv")


################################################################################################
#            calculate mean and standard deviation of the three replicates                     #
################################################################################################

Uptake.ratio.np$ratio_long<- apply(Uptake.ratio.np[,2:4], 1, mean)

Uptake.ratio.np$ratio_short<- apply(Uptake.ratio.np[,5:7], 1, mean)

Uptake.ratio.np$sd_long<- apply(Uptake.ratio.np[,2:4], 1, sd)

Uptake.ratio.np$sd_short<- apply(Uptake.ratio.np[,5:7], 1, sd)

m1<- apply(Uptake.ratio.np[,2:9], MARGIN = 2, mean)

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

write.csv(Uptake.ratio.np, "./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")


##################################################################################################
# make a normalized depth dataframe and a pseudocount uptake ratio                               #
##################################################################################################

# make a normalized depth dataframe
norm.ratios.np<-  data.table(pos = remix13.UP01$pos, UP01 = remix13.UP01$UP01.n, 
                             UP02 = remix13.UP02$UP02.n, UP03 = remix13.UP03$UP03.n,
                             UP07 = remix15.UP07$UP07.n, UP08 = remix15.UP08$UP08.n, 
                             UP09 = remix15.UP09$UP09.n, UP13.UP01 = remix13.UP01$input.n, 
                             UP13.UP02 = remix13.UP02$input.n, UP13.UP03 = remix13.UP03$input.n, 
                             UP15.UP07 = remix15.UP07$input.n, UP15.UP08 = remix15.UP08$input.n, 
                             UP15.UP09 = remix15.UP09$input.n) 
# save
write.csv(norm.ratios.np, "./datasets/final_datasets/norm.ratios.np.corrected.csv")



