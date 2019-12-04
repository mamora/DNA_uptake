################################################################################################
#     This script takes the raw bed files from the sequence handling scripts and calculates    #
#     normalized uptake ratios and well as normalized depth files for GG                       #
################################################################################################
#     In this case I will use bed files from samples aligned to concatenated Rd/GG and Rd/GG   #
#      and where low quality reads (mapping quality  = 0 ) have been removed                   #
################################################################################################

# load package
library(data.table)


#######################################################
# load input and donor bed files and calculate uptake #
#######################################################

UP04 <- read.delim("./raw_depth/map_q/both_genomes/UP04.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP04) <- c("genome","pos","UP04") # add headers

UP05 <- read.delim("./raw_depth/map_q/both_genomes/UP05.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP05) <- c("genome","pos","UP05") # add headers

UP06 <- read.delim("./raw_depth/map_q/both_genomes/UP06.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP06) <- c("genome","pos","UP06") # add headers

UP10 <- read.delim("./raw_depth/map_q/both_genomes/UP10.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP10) <- c("genome","pos","UP10") # add headers

UP11 <- read.delim("./raw_depth/map_q/both_genomes/UP11.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP11) <- c("genome","pos","UP11") # add headers

UP12 <- read.delim("./raw_depth/map_q/both_genomes/UP12.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP12) <- c("genome","pos","UP12") # add headers

UP14 <- read.delim("./raw_depth/map_q/both_genomes/UP14.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP14) <- c("genome","pos","UP14") # add headers

UP16<- read.delim("./raw_depth/map_q/both_genomes/UP16.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP16)<- c("genome","pos","UP16") # add headers

UP14.new <- read.delim("./raw_depth/map_q/both_genomes/UP14new.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP14.new) <- c("genome","pos","UP14.new") # add headers

UP16.new<- read.delim("./raw_depth/map_q/both_genomes/UP16new.rdgg_19_map_depth.bed", header=FALSE)  
colnames(UP16.new)<- c("genome","pos","UP16.new") # add headers



##########################################################################################
#     generate a table of raw-depth from all reads of all samples aligned to GG only     #
##########################################################################################


#get only depth of reads aligned to GG
GG_only_raw_reads <- UP04[grep("L_contig000001", UP04$genome), 1:dim(UP04)[2]] 

# remove 1 column
GG_only_raw_reads<- GG_only_raw_reads[,2:3]

# add the rest of samples
GG_only_raw_reads$UP05<- UP05[grep("L_contig000001", UP05$genome), "UP05"] 
GG_only_raw_reads$UP06<- UP06[grep("L_contig000001", UP06$genome), "UP06"] 
GG_only_raw_reads$UP10<- UP10[grep("L_contig000001", UP10$genome), "UP10"] 
GG_only_raw_reads$UP11<- UP11[grep("L_contig000001", UP11$genome), "UP11"] 
GG_only_raw_reads$UP12<- UP12[grep("L_contig000001", UP12$genome), "UP12"] 
GG_only_raw_reads$UP14<- UP14[grep("L_contig000001", UP14$genome), "UP14"] 
GG_only_raw_reads$UP16<- UP16[grep("L_contig000001", UP16$genome), "UP16"] 
GG_only_raw_reads$UP14.new<- UP14.new[grep("L_contig000001", UP14.new$genome), "UP14.new"] 
GG_only_raw_reads$UP16.new<- UP16.new[grep("L_contig000001", UP16.new$genome), "UP16.new"] 

write.csv(GG_only_raw_reads, "raw_depth/datasets/map_q/GG_only_raw_reads_map_contam_correct.csv")

##########################################################################################
#     generate a table of raw-depth from all reads of all samples aligned to Rd only     #
##########################################################################################


#get only depth of reads aligned to Rd
Rd_only_raw_reads <- UP04[grep("NC_000907.1", UP04$genome), 1:dim(UP04)[2]] 

# remove 1 column
Rd_only_raw_reads <- Rd_only_raw_reads[,2:3]

# add the rest of samples
Rd_only_raw_reads$UP05 <- UP05[grep("NC_000907.1", UP05$genome), "UP05"] 
Rd_only_raw_reads$UP06 <- UP06[grep("NC_000907.1", UP06$genome), "UP06"] 
Rd_only_raw_reads$UP10 <- UP10[grep("NC_000907.1", UP10$genome), "UP10"] 
Rd_only_raw_reads$UP11 <- UP11[grep("NC_000907.1", UP11$genome), "UP11"] 
Rd_only_raw_reads$UP12 <- UP12[grep("NC_000907.1", UP12$genome), "UP12"] 
Rd_only_raw_reads$UP14 <- UP14[grep("NC_000907.1", UP14$genome), "UP14"] 
Rd_only_raw_reads$UP16 <- UP16[grep("NC_000907.1", UP16$genome), "UP16"] 
Rd_only_raw_reads$UP14.new<- UP14.new[grep("NC_000907.1", UP14.new$genome), "UP14.new"] 
Rd_only_raw_reads$UP16.new<- UP16.new[grep("NC_000907.1", UP16.new$genome), "UP16.new"] 


write.csv(Rd_only_raw_reads, "raw_depth/datasets/map_q/Rd_onlyGG_raw_reads_map_contam_correct.csv")


##########################################################################################
#                      Normalize to depth to depth/average depth per position            #
##########################################################################################

norm_s <- function(colm = GG_only_raw_reads$UP04) {
  norm <- colm/(sum(colm)/length(colm))
  return(norm)
}

table1 <- data.table(pos = GG_only_raw_reads$pos)

for (i in 2:length(GG_only_raw_reads)){
  n <- norm_s(GG_only_raw_reads[,i])
  table1 <- cbind(table1, n)
}

name <- colnames(GG_only_raw_reads)

colnames(table1) <- name 

write.csv(table1, "raw_depth/datasets/map_q/GG_only_norm_reads_map_contam_correct.csv")

summary(table1)



table2 <- data.table(pos = Rd_only_raw_reads$pos)

for (i in 2: length(Rd_only_raw_reads)){
  n <- norm_s(Rd_only_raw_reads[,i])
  table2 <- cbind(table2, n)
}

name <- colnames(Rd_only_raw_reads)

colnames(table2) <- name 

write.csv(table2, "raw_depth/datasets/map_q/Rd_onlyGG_norm_reads_map_contam_correct.csv")


##########################################################################################
#                             Calculate Uptake ratio for GG                              #
##########################################################################################


GG_norm_reads_contam_correct <- fread("raw_depth/datasets/map_q/GG_only_norm_reads_map_contam_correct.csv")

GG_norm_reads_contam_correct$V1 <- NULL

# divide donor normalize depth by input normalze depth
u4 <- GG_norm_reads_contam_correct$UP04/GG_norm_reads_contam_correct$UP14.new
u5 <- GG_norm_reads_contam_correct$UP05/GG_norm_reads_contam_correct$UP14.new
u6 <- GG_norm_reads_contam_correct$UP06/GG_norm_reads_contam_correct$UP14.new
u10 <- GG_norm_reads_contam_correct$UP10/GG_norm_reads_contam_correct$UP16.new
u11 <- GG_norm_reads_contam_correct$UP11/GG_norm_reads_contam_correct$UP16.new
u12 <- GG_norm_reads_contam_correct$UP12/GG_norm_reads_contam_correct$UP16.new

GG_contam_correct_uptake_ratio <- data.table(pos = GG_norm_reads_contam_correct$pos, UP04 = u4,
                                             UP05 = u5, UP06 = u6, UP10 = u10, UP11 = u11, UP12 = u12)

# calculate mean and standard deviation
GG_contam_correct_uptake_ratio$ratio_large <- apply(GG_contam_correct_uptake_ratio[,2:4], 1, mean)

GG_contam_correct_uptake_ratio$ratio_short <- apply(GG_contam_correct_uptake_ratio[,5:7], 1, mean)

GG_contam_correct_uptake_ratio$sd_large <- apply(GG_contam_correct_uptake_ratio[,2:4], 1, sd)

GG_contam_correct_uptake_ratio$sd_short <- apply(GG_contam_correct_uptake_ratio[,5:7], 1, sd)


write.csv(GG_contam_correct_uptake_ratio, "raw_depth/datasets/map_q/GG_contam_correct_map_uptake_ratio.csv")


##############################################
#   add flags for input low coverage areas   #
##############################################


GG_only_raw_reads <- fread("raw_depth/datasets/map_q/GG_only_raw_reads_map_contam_correct.csv")
GG_contam_correct_uptake_ratio <- fread("raw_depth/datasets/map_q/GG_contam_correct_map_uptake_ratio.csv")



length(which(GG_only_raw_reads$UP16.new == 0))/length((GG_only_raw_reads$UP16.new))

flag_short <- GG_only_raw_reads$UP16.new >= 20  

flag_large <- GG_only_raw_reads$UP14.new >= 20   

GG_contam_correct_uptake_ratio$flag_short <- flag_short
GG_contam_correct_uptake_ratio$flag_large <- flag_large

write.csv(GG_contam_correct_uptake_ratio, "raw_depth/datasets/map_q/GG_contam_correct_map_uptake_ratio.csv")


plot(GG_contam_correct_uptake_ratio$pos[1:50000], GG_contam_correct_uptake_ratio$ratio_short[1:50000])



##############################################
#             table of coverages             #
##############################################


cut1<- which(GG_contam_correct_uptake_ratio$smooth_large_ratio < 0.25)
round(((length(cut1)/length(GG_contam_correct_uptake_ratio$pos)) * 100), digits = 2)

test<- GG_only_raw_reads$UP14.new[cut1]


round(((length(which(test <= 10))/length(test)) * 100), digits = 3)
round(((length(which(test > 10 & test <= 20 ))/length(test)) * 100), digits = 3)
round(((length(which(test > 20 & test <= 50))/length(test)) * 100), digits = 3)
round(((length(which(test > 50 & test <= 100))/length(test)) * 100), digits = 3)
round(((length(which(test > 100))/length(test)) * 100), digits = 3)




cut1<- which(GG_contam_correct_uptake_ratio$smooth_small_ratio > 5)
round(((length(cut1)/length(GG_contam_correct_uptake_ratio$pos)) * 100), digits = 2)

test<- GG_only_raw_reads$UP16.new[cut1]


round(((length(which(test <= 10))/length(test)) * 100), digits = 3)
round(((length(which(test > 10 & test <= 20 ))/length(test)) * 100), digits = 3)
round(((length(which(test > 20 & test <= 50))/length(test)) * 100), digits = 3)
round(((length(which(test > 50 & test <= 100))/length(test)) * 100), digits = 3)
round(((length(which(test > 100))/length(test)) * 100), digits = 3)


