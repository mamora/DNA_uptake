################################################################################################
#     This script takes the raw bed files from the sequence handling scripts and calculates    #
#     normalized uptake ratios and well as normalized depth files for NP and GG genome         #
################################################################################################
#     In this case I will use bed files from samples aligned to concatenated Rd/NP and Rd/GG   #
#      and where low quality reads (mapping quality  = 0 ) have been removed                   #
################################################################################################

library(data.table)


#load UP01 depth  
# load input and donor bed files and calculate uptake for UP01
UP01 <- read.delim("./raw_depth/map_q/both_genomes/UP01.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP01) <- c("genome","pos","UP01") # add headers

UP02 <- read.delim("./raw_depth/map_q/both_genomes/UP02.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP02) <- c("genome","pos","UP02") # add headers

UP03 <- read.delim("./raw_depth/map_q/both_genomes/UP03.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP03) <- c("genome","pos","UP03") # add headers

UP07 <- read.delim("./raw_depth/map_q/both_genomes/UP07.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP07) <- c("genome","pos","UP07") # add headers

UP08 <- read.delim("./raw_depth/map_q/both_genomes/UP08.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP08) <- c("genome","pos","UP08") # add headers

UP09 <- read.delim("./raw_depth/map_q/both_genomes/UP09.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP09) <- c("genome","pos","UP09") # add headers

UP13 <- read.delim("./raw_depth/map_q/both_genomes/UP13.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP13) <- c("genome","pos","UP13") # add headers

UP15<- read.delim("./raw_depth/map_q/both_genomes/UP15.rdnp_19_map_depth.bed", header=FALSE)  
colnames(UP15)<- c("genome","pos","UP15") # add headers


##########################################################################################
#     generate a table of raw-depth from all reads of all samples aligned to NP only     #
##########################################################################################


#get only depth of reads aligned to NP
NP_only_raw_reads <- UP01[grep("NC_007146.2", UP01$genome), 1:dim(UP01)[2]] 

# remove 1 column
NP_only_raw_reads<- NP_only_raw_reads[,2:3]

# add the rest of samples
NP_only_raw_reads$UP02<- UP02[grep("NC_007146.2", UP02$genome), "UP02"] 
NP_only_raw_reads$UP03<- UP03[grep("NC_007146.2", UP03$genome), "UP03"] 
NP_only_raw_reads$UP07<- UP07[grep("NC_007146.2", UP07$genome), "UP07"] 
NP_only_raw_reads$UP08<- UP08[grep("NC_007146.2", UP08$genome), "UP08"] 
NP_only_raw_reads$UP09<- UP09[grep("NC_007146.2", UP09$genome), "UP09"] 
NP_only_raw_reads$UP13<- UP13[grep("NC_007146.2", UP13$genome), "UP13"] 
NP_only_raw_reads$UP15<- UP15[grep("NC_007146.2", UP15$genome), "UP15"] 

write.csv(NP_only_raw_reads, "raw_depth/datasets/map_q/NP_only_raw_reads_map_contam_correct.csv")

##########################################################################################
#     generate a table of raw-depth from all reads of all samples aligned to Rd only     #
##########################################################################################


#get only depth of reads aligned to Rd
Rd_only_raw_reads <- UP01[grep("NC_000907.1", UP01$genome), 1:dim(UP01)[2]] 

# remove 1 column
Rd_only_raw_reads <- Rd_only_raw_reads[,2:3]

# add the rest of samples
Rd_only_raw_reads$UP02 <- UP02[grep("NC_000907.1", UP02$genome), "UP02"] 
Rd_only_raw_reads$UP03 <- UP03[grep("NC_000907.1", UP03$genome), "UP03"] 
Rd_only_raw_reads$UP07 <- UP07[grep("NC_000907.1", UP07$genome), "UP07"] 
Rd_only_raw_reads$UP08 <- UP08[grep("NC_000907.1", UP08$genome), "UP08"] 
Rd_only_raw_reads$UP09 <- UP09[grep("NC_000907.1", UP09$genome), "UP09"] 
Rd_only_raw_reads$UP13 <- UP13[grep("NC_000907.1", UP13$genome), "UP13"] 
Rd_only_raw_reads$UP15 <- UP15[grep("NC_000907.1", UP15$genome), "UP15"] 


write.csv(Rd_only_raw_reads, "raw_depth/datasets/map_q/Rd_only_raw_reads_map_contam_correct.csv")


##########################################################################################
#                      Normalize to depth to depth/average depth per position            #
##########################################################################################

norm_s <- function(colm = NP_only_raw_reads$UP01) {
  norm <- colm/(sum(colm)/length(colm))
  return(norm)
}

table1 <- data.table(pos = NP_only_raw_reads$pos)

for (i in 2:length(NP_only_raw_reads)){
  n <- norm_s(NP_only_raw_reads[,i])
  table1 <- cbind(table1, n)
}

name <- colnames(NP_only_raw_reads)

colnames(table1) <- name 

write.csv(table1, "raw_depth/datasets/map_q/NP_only_norm_reads_map_contam_correct.csv")

summary(table1)



table2 <- data.table(pos = Rd_only_raw_reads$pos)

for (i in 2: length(Rd_only_raw_reads)){
  n <- norm_s(Rd_only_raw_reads[,i])
  table2 <- cbind(table2, n)
}

name <- colnames(Rd_only_raw_reads)

colnames(table2) <- name 

write.csv(table2, "raw_depth/datasets/map_q/Rd_only_norm_reads_map_contam_correct.csv")


##########################################################################################
#                             Calculate Uptake ratio for NP                              #
##########################################################################################


NP_norm_reads_contam_correct <- fread("raw_depth/datasets/map_q/NP_only_norm_reads_map_contam_correct.csv")

NP_norm_reads_contam_correct$V1 <- NULL

# divide donor normalize depth by input normalze depth
u1 <- NP_norm_reads_contam_correct$UP01/NP_norm_reads_contam_correct$UP13
u2 <- NP_norm_reads_contam_correct$UP02/NP_norm_reads_contam_correct$UP13
u3 <- NP_norm_reads_contam_correct$UP03/NP_norm_reads_contam_correct$UP13
u7 <- NP_norm_reads_contam_correct$UP07/NP_norm_reads_contam_correct$UP15
u8 <- NP_norm_reads_contam_correct$UP08/NP_norm_reads_contam_correct$UP15
u9 <- NP_norm_reads_contam_correct$UP09/NP_norm_reads_contam_correct$UP15

NP_contam_correct_uptake_ratio <- data.table(pos = NP_norm_reads_contam_correct$pos, UP01 = u1,
                                             UP02 = u2, UP03 = u3, UP07 = u7, UP08 = u8, UP09 = u9)

# calculate mean and standard deviation
NP_contam_correct_uptake_ratio$ratio_large <- apply(NP_contam_correct_uptake_ratio[,2:4], 1, mean)

NP_contam_correct_uptake_ratio$ratio_short <- apply(NP_contam_correct_uptake_ratio[,5:7], 1, mean)

NP_contam_correct_uptake_ratio$sd_large <- apply(NP_contam_correct_uptake_ratio[,2:4], 1, sd)

NP_contam_correct_uptake_ratio$sd_short <- apply(NP_contam_correct_uptake_ratio[,5:7], 1, sd)


write.csv(NP_contam_correct_uptake_ratio, "raw_depth/datasets/map_q/NP_contam_correct_map_uptake_ratio.csv")


##############################################
#   add flags for input low coverage areas   #
##############################################

NP_only_raw_reads <- fread("./raw_depth/datasets/map_q/NP_only_raw_reads_map_contam_correct.csv")
NP_contam_correct_uptake_ratio <- fread("raw_depth/datasets/map_q/NP_contam_correct_map_uptake_ratio.csv")



flag_short <- NP_only_raw_reads$UP15 > 20  

flag_large <- NP_only_raw_reads$UP13 > 20   

NP_contam_correct_uptake_ratio$flag_short <- flag_short
NP_contam_correct_uptake_ratio$flag_large <- flag_large

write.csv(NP_contam_correct_uptake_ratio, "raw_depth/datasets/map_q/NP_contam_correct_map_uptake_ratio.csv")


plot(NP_contam_correct_uptake_ratio$pos[10000:20000], NP_contam_correct_uptake_ratio$ratio_short[10000:20000])


summary(NP_only_raw_reads$UP15)
View(NP_only_raw_reads[1230384:1231349,])

summary(NP_only_raw_reads$UP15)

hist(NP_only_raw_reads$UP15, breaks=seq(0,2000,by=10), xlab = "raw reads coverage", freq=FALSE, axes = FALSE, main = "histogram of small fragment input raw reads density")
axis(side = 1, at = c(seq(0,2000,by=100)))
axis(side = 2, at = c(seq(0,0.004,by=0.001)))



##############################################
#             table of coverages             #
##############################################

cut0<- which(NP_only_raw_reads$UP13 == 0)

round(((length(cut0)/length(NP_contam_correct_uptake_ratio$pos)) * 100), digits = 2)

test<- NP_contam_correct_uptake_ratio[cut0,] 





cut1<- which(NP_contam_correct_uptake_ratio$smooth_large_ratio >= 2)
round(((length(cut1)/length(NP_contam_correct_uptake_ratio$pos)) * 100), digits = 2)

test<- NP_only_raw_reads$UP13[cut1]


round(((length(which(test <= 10))/length(test)) * 100), digits = 2)
round(((length(which(test > 10 & test <= 20 ))/length(test)) * 100), digits = 2)
round(((length(which(test > 20 & test <= 50))/length(test)) * 100), digits = 2)
round(((length(which(test > 50 & test <= 100))/length(test)) * 100), digits = 2)
round(((length(which(test > 100))/length(test)) * 100), digits = 2)





