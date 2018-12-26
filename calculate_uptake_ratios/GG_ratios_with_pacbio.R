################################################################################################
#     This script takes the raw bed files from the sequence handling scripts and calculates    #
#     normalized uptake ratios and well as normalized depth files for GG genome                #
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

###################################     UP04       ######################################

# load input and donor bed files and calculate uptake for UP04
remix14.UP04<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/REMIX14.UP04.gg.pacbio.depth.bed") #load new input samples fro UP01   
colnames(remix14.UP04)<- c("genome","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up04<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/UP04.pilon_pittgg_pacBio.depth.bed")
colnames(up04)<- c("genome","pos","depth") # add headers
remix14.UP04$UP04<- up04$depth
remix14.UP04$input.mean<- (remix14.UP04$depth.1+remix14.UP04$depth.2+remix14.UP04$depth.3+remix14.UP04$depth.4)/4

# Normalize the raw reads to depth per millon reads
up14.up04.mapped<-  sum(remix14.UP04$input.mean)/length(remix14.UP04$input.mean)  
up04.all.mapped<- sum(remix14.UP04$UP04)/length(remix14.UP04$input.mean)
remix14.UP04$UP04.n<- (remix14.UP04$UP04)/up04.all.mapped
remix14.UP04$input.n<- (remix14.UP04$input.mean)/up14.up04.mapped
remix14.UP04$UP04.ratio<- remix14.UP04$UP04.n/remix14.UP04$input.n

###################################     UP05       ######################################

# load input and donor bed files and calculate uptake for UP05
remix14.UP05<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/REMIX14.UP05.gg.pacbio.depth.bed") #load new input samples fro UP01   
colnames(remix14.UP05)<- c("genome","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up05<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/UP05.pilon_pittgg_pacBio.depth.bed")
colnames(up05)<- c("genome","pos","depth") # add headers
remix14.UP05$UP05<- up05$depth
remix14.UP05$input.mean<- (remix14.UP05$depth.1+remix14.UP05$depth.2+remix14.UP05$depth.3+remix14.UP05$depth.4)/4


# Normalize the raw reads to depth per millon reads
up14.up05.mapped<- sum(remix14.UP05$input.mean)/length(remix14.UP04$input.mean)   # total mapped reads, from sambamba flagstats files 
up05.all.mapped<- sum(remix14.UP05$UP05)/length(remix14.UP04$input.mean)   # total mapped reads from summary table
remix14.UP05$UP05.n<- (remix14.UP05$UP05)/up05.all.mapped
remix14.UP05$input.n<- (remix14.UP05$input.mean)/up14.up05.mapped
remix14.UP05$UP05.ratio<- remix14.UP05$UP05.n/remix14.UP05$input.n

###################################     UP06       ######################################

# load input and donor bed files and calculate uptake for UP06
remix14.UP06<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/REMIX14.UP06.gg.pacbio.depth.bed") #load new input samples fro UP01   
colnames(remix14.UP06)<- c("genome","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up06<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/UP06.pilon_pittgg_pacBio.depth.bed")
colnames(up06)<- c("genome","pos","depth") # add headers
remix14.UP06$UP06<- up06$depth
remix14.UP06$input.mean<- (remix14.UP06$depth.1+remix14.UP06$depth.2+remix14.UP06$depth.3+remix14.UP06$depth.4)/4
#########################################################################

# Normalize the raw reads to depth per millon reads
up14.up06.mapped<- sum(remix14.UP06$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up06.all.mapped<- sum(remix14.UP06$UP06)/length(remix14.UP06$input.mean)  # total mapped reads from summary table
remix14.UP06$UP06.n<- (remix14.UP06$UP06)/up06.all.mapped
remix14.UP06$input.n<- (remix14.UP06$input.mean)/up14.up06.mapped
remix14.UP06$UP06.ratio<- remix14.UP06$UP06.n/remix14.UP06$input.n

# save in dataframes
raw.depth.samples.gg<- data.table(pos = remix14.UP04$pos, UP04 = remix14.UP04$UP04, UP05 = remix14.UP05$UP05, UP06 = remix14.UP06$UP06) 
raw.depth.inputs.gg<-  data.table(pos = remix14.UP04$pos, UP14.UP04 = remix14.UP04$input.mean, UP14.UP05 = remix14.UP05$input.mean, UP14.UP06 = remix14.UP06$input.mean) 
Uptake.ratio.gg<-  data.table(pos = remix14.UP04$pos, UP04.ratio = remix14.UP04$UP04.ratio, UP05.ratio = remix14.UP05$UP05.ratio, UP06.ratio = remix14.UP06$UP06.ratio) 

#################################      UP10       ############################################

# load input and donor bed files and calculate uptake for UP10
remix16.UP10<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/REMIX16.UP10.gg.pacbio.depth.bed") #load new input samples fro UP01   
colnames(remix16.UP10)<- c("genome","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up10<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/UP10.pilon_pittgg_pacBio.depth.bed")
colnames(up10)<- c("genome","pos","depth") # add headers
remix16.UP10$UP10<- up10$depth
remix16.UP10$input.mean<- (remix16.UP10$depth.1+remix16.UP10$depth.2+remix16.UP10$depth.3+remix16.UP10$depth.4)/4

# Normalize the raw reads to depth per millon reads
up16.up10.mapped<- sum(remix16.UP10$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up10.all.mapped<- sum(remix16.UP10$UP10)/length(remix14.UP06$input.mean) # total mapped reads from summary table
remix16.UP10$UP10.n<- (remix16.UP10$UP10)/up10.all.mapped
remix16.UP10$input.n<- (remix16.UP10$input.mean)/up16.up10.mapped
remix16.UP10$UP10.ratio<- remix16.UP10$UP10.n/remix16.UP10$input.n

# save in dataframes
raw.depth.samples.gg$UP10<- remix16.UP10$UP10
raw.depth.inputs.gg$UP16.UP10<- remix16.UP10$input.mean
Uptake.ratio.gg$UP10.ratio<- remix16.UP10$UP10.ratio

######################################      UP11        ###################################################

# load input and donor bed files and calculate uptake for UP11
remix16.UP11<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/REMIX16.UP11.gg.pacbio.depth.bed") #load new input samples fro UP01   
colnames(remix16.UP11)<- c("genome","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up11<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/UP11.pilon_pittgg_pacBio.depth.bed")
colnames(up11)<- c("genome","pos","depth") # add headers
remix16.UP11$UP11<- up11$depth
remix16.UP11$input.mean<- (remix16.UP11$depth.1+remix16.UP11$depth.2+remix16.UP11$depth.3+remix16.UP11$depth.4)/4

# Normalize the raw reads to depth per millon reads
up16.up11.mapped<- sum(remix16.UP11$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up11.all.mapped<- sum(remix16.UP11$UP11)/length(remix14.UP06$input.mean)  # total mapped reads from summary table
remix16.UP11$UP11.n<- (remix16.UP11$UP11)/up11.all.mapped
remix16.UP11$input.n<- (remix16.UP11$input.mean)/up16.up11.mapped
remix16.UP11$UP11.ratio<- remix16.UP11$UP11.n/remix16.UP11$input.n

raw.depth.samples.gg$UP11<- remix16.UP11$UP11
raw.depth.inputs.gg$UP16.UP11<- remix16.UP11$input.mean
Uptake.ratio.gg$UP11.ratio<- remix16.UP11$UP11.ratio

######################################      UP12        ############################################

# load input and donor bed files and calculate uptake for UP12
remix16.UP12<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/REMIX16.UP12.gg.pacbio.depth.bed") #load new input samples fro UP01   
colnames(remix16.UP12)<- c("genome","pos","depth.1","depth.2","depth.3","depth.4") # add headers
up12<- fread("C:/Users/marcelo/Dropbox/uptake/PittGG_new_input_raw_reads/depth/UP12.pilon_pittgg_pacBio.depth.bed")
colnames(up12)<- c("genome","pos","depth") # add headers
remix16.UP12$UP12<- up12$depth
remix16.UP12$input.mean<- (remix16.UP12$depth.1+remix16.UP12$depth.2+remix16.UP12$depth.3+remix16.UP12$depth.4)/4

# Normalize the raw reads to depth per millon reads
up16.up12.mapped<- sum(remix16.UP12$input.mean)/length(remix14.UP06$input.mean)  # total mapped reads, from sambamba flagstats files 
up12.all.mapped<- sum(remix16.UP12$UP12)/length(remix14.UP06$input.mean) # total mapped reads from summary table
remix16.UP12$UP12.n<- (remix16.UP12$UP12)/up12.all.mapped
remix16.UP12$input.n<- (remix16.UP12$input.mean)/up16.up12.mapped
remix16.UP12$UP12.ratio<- remix16.UP12$UP12.n/remix16.UP12$input.n

raw.depth.samples.gg$UP12<- remix16.UP12$UP12
raw.depth.inputs.gg$UP16.UP12<- remix16.UP12$input.mean
Uptake.ratio.gg$UP12.ratio<- remix16.UP12$UP12.ratio


write.csv(raw.depth.samples.gg, "./datasets/final_datasets/pacbio.raw.depth.samples.gg.csv")
write.csv(raw.depth.inputs.gg, "./datasets/final_datasets/pacbio.raw.depth.inputs.gg.csv")

################################################################################################
#            calculate mean and standard deviation of the three replicates                     #
################################################################################################

Uptake.ratio.gg$ratio_long<- apply(Uptake.ratio.gg[,2:4], 1, mean)

Uptake.ratio.gg$ratio_short<- apply(Uptake.ratio.gg[,5:7], 1, mean)

Uptake.ratio.gg$sd_long<- apply(Uptake.ratio.gg[,2:4], 1, sd)

Uptake.ratio.gg$sd_short<- apply(Uptake.ratio.gg[,5:7], 1, sd)

# function for normalizing predicted uptake
norm<- function (data = data){
  h_pe<- mean(data, na.rm = TRUE)
  s_pe<- (data * 1)/h_pe
  return(s_pe)
}

# normalize predicted uptake to a mean of 1
Uptake.ratio.gg$ratio_short<- norm(data = Uptake.ratio.gg$ratio_short)
Uptake.ratio.gg$ratio_long<- norm(data = Uptake.ratio.gg$ratio_long)

mean(Uptake.ratio.gg$ratio_short, na.rm = TRUE)
mean(Uptake.ratio.gg$ratio_long, na.rm = TRUE)

write.csv(Uptake.ratio.gg, "./datasets/final_datasets/DNA_uptake/Uptake.ratio.gg.csv")

##################################################################################################
# make a normalized depth dataframe and a pseudocount uptake ratio                               #
##################################################################################################

# make a normalized depth dataframe
norm.ratios.gg<-  data.table(pos = remix14.UP04$pos, UP04 = remix14.UP04$UP04.n, 
                             UP05 = remix14.UP05$UP05.n, UP06 = remix14.UP06$UP06.n,
                             UP10 = remix16.UP10$UP10.n, UP11 = remix16.UP11$UP11.n,
                             UP12 = remix16.UP12$UP12.n, UP14.UP04 = remix14.UP04$input.n, 
                             UP14.UP05 = remix14.UP05$input.n, UP14.UP06 = remix14.UP06$input.n,
                             UP16.UP10 = remix16.UP10$input.n, UP16.UP11 = remix16.UP11$input.n, 
                             UP16.UP12 = remix16.UP12$input.n) 

# save
write.csv(norm.ratios.gg, "./datasets/final_datasets/norm.ratios.gg.corrected.csv")
