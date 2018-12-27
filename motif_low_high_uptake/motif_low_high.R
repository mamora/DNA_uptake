####################################################################################
#     This script will build motifs for USS with low/high uptake and scores        #
####################################################################################


# load dataframes
library(here)
library(plyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)

folder.name <- "./Marcelo_paper_figures_new_order/supplementary/suppl_9/"  
#read uptake file
Uptake.ratio.np <- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")      
# load USS scores
np.USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")
# list of USS10
Up.USS.np.10.list <- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) 
Up.USS.np.9.5.list<- fread(here::here("datasets/final_datasets","Up.USS.np.9.5.list.csv")) #save file

# load functions
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")


##################################################################################################
#                     make a dataframe of USS that are far away from the closest USS by 800 bp   #
##################################################################################################

#calculate closest distance to closest USS for each position
close<- sapply(Up.USS.np.9.5.list$keypos,dist.USS, USS.genome.c = Up.USS.np.10.list$keypos) 

far.uss.9.5<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$next.uss >= 800),]


##################################################################################################
#                     calculate uptake of 31 bases from focal position                           # 
##################################################################################################

# make a dataframe uptake ratios and scores with first positions at the end to circularize the genome
data.np<- c(Uptake.ratio.np$ratio_short, Uptake.ratio.np$ratio_short[1:30]) 

# calculate mean ratio of 31 positions left of each pos
mean.ratio<- rollapply(data = data.np, width = 31, align = "left", by = 1, FUN = mean) 
joined<- data.frame(pos = Uptake.ratio.np$pos, mean.ratio = mean.ratio, max.score = np.USS.scores$max)
joined$strand<- np.USS.scores$strand # add strand of highest score per focal position

################################# extract positions with uss scores aroung 10 and uptake ratios high and low   ####################################

##################################################################################################
#                     subset USS by uptake ratio and extract DNA sequence of those USS           # 
##################################################################################################


# filter uss by USS score
Up.USS.np.9.5.list$mean.ratio<- joined$mean.ratio[Up.USS.np.9.5.list$USS.pos]

high.ratio<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$mean.ratio >=  3),]
mid.ratio<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$mean.ratio >= 0.5 & Up.USS.np.9.5.list$mean.ratio < 3),]
low.ratio<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$mean.ratio < 0.5),]


# get NP genomic sequence
# Name the sites file, and genome reference file
fasta.np    <- "datasets/sequences/np.pilon.fasta" # Reference fasta to the NP genome
genome.np <- read.genome(fasta.np) # read fasta with seqinr, then massage

# function to extract sequence of USS
make.seq.list<- function(low.uss, save){
  for (i in 1:length(low.uss$USS.pos)){
    if (low.uss$strand[i] == "w") { 
      site = low.uss$USS.pos[i]
      w<- as.character(genome.np[[1]][(site):(site + 30)])
      write.fasta(sequences = w, names = low.uss$USS.pos[i], open = "a", file.out = save)
      low.uss.10<- c(low.uss.10,paste(w, collapse = ""))
    } else {
      site = low.uss$USS.pos[i]
      g<- genome.np[[1]][(site + 30):(site)]
      c<-  as.character(revalue(g, c("a"="t", "t"="a","g"="c","c"="g","n" = "n")))
      write.fasta(sequences = c, names = low.uss$USS.pos[i], open = "a", file.out = save)
      low.uss.10<- c(low.uss.10,paste(c, collapse = ""))
    }  
  }  
  return(low.uss.10)
}


# extract sequences of USSs subset above
low.uss.10<- NULL
save<- paste(folder.name,"motif_high","more_3.fa", sep = "")
uss.h.seq<- make.seq.list(high.ratio, save)

low.uss.10<- NULL
save<- paste(folder.name,"motif_mid","more_0_5_less_3.fa", sep = "")
uss.m.seq<- make.seq.list(mid.ratio, save)

low.uss.10<- NULL
save<- paste(folder.name,"motif_low","less_0_5.fa", sep = "")
uss.l.seq<- make.seq.list(low.ratio, save)


##################################################################################################
#                     subset USS by uptake ratio and extract DNA sequence of those USS           # 
##################################################################################################

# Name the sites file, and genome reference file
sites<- uss.h.seq  
sites1<- uss.m.seq
sites2<- uss.l.seq

# download seqLogo package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("seqLogo", version = "3.8")

# load packages for building motifs
library(seqLogo)
library(RColorBrewer)

# load functions to build PWM matrix, functions in this file are a modified version of seqlogo 
# with a customized y-axis  
source(here::here("helper_functions","seqLogoFunction_y_axis.R"))

# function to make a matrix out of every base in the file
sites2matrix_alt<- function(file="Hin2206sites.txt"){
  sites    <- file                         # rename file
  sites    <- tolower(sites)               # make everything lowercase
  nsites   <- length(sites)                # how many sites?
  nbases   <- nchar(sites[1])              # how long is the first site?
  if(all(nchar(sites) == nbases) == T){    # return matrix, if sites all same length
    sitestr <- stringsplit(sites)        # vectorize
    goodchars <- which(sitestr %in% chars[-length(chars)]) # index ACGT pos, 
    sitestr[-goodchars] <- chars[length(chars)]            # Nify non-ACGT
    sitemat <- matrix(sitestr, nrow=nsites, ncol=nbases, byrow=TRUE) # matrix
    return(sitemat)   # output results
  }else{return(NA)} # return NA, if sites not all same length
}


# function to make a logo  
make.logo<- function(sites6, save){
  sitemat <- sites2matrix_alt(sites6)    # set of sites
  pcm     <- matrix2pcm(sitemat)    # count matrix
  pfm     <- pcm2pfm(pcm, pseudo=0) # frequency matrix
  #make the PPM or PFM matrix used for the logo
  p<- makePWM(pfm)
  #save
  pdf(file= save, width = 12, height = 7)
  seqLogo(p)
  dev.off()
}

# build logos
make.logo(sites, save = paste(folder.name,"motif_high.test2","more_3.pdf", sep = ""))
make.logo(sites1, save = paste(folder.name,"motif_mid","more_0_5_less_3.pdf", sep = ""))
make.logo(sites2, save = paste(folder.name,"motif_low","less_0_5.pdf", sep = ""))



#################################################################################
#                    Extract outlier USS sequences                              #
#################################################################################

# get stats from outliers taken from table S4 of the manuscript
outlier_higher<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$USS.pos %in% c(547194,532339, 977257, 
                                                                            1593554, 1788333)),]

outlier_lower<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$keypos %in% c(508994,552251, 587402, 
                                                                           701942, 836753,976517,
                                                                           1061432,112312,1630026,
                                                                           1653648,1828950))]


low.uss.10<- NULL
save<- paste(folder.name,"outlier_higher","seq.fa", sep = "")
outlier_higher.seq<- make.seq.list(outlier_higher, save)

low.uss.10<- NULL
save<- paste(folder.name,"outlier_lower","seq.fa", sep = "")
outlier_lower.seq<- make.seq.list(outlier_lower, save)

pos<- c(547194,532339, 977257,1593554, 1788333)
seq<- outlier_higher.seq

# get sequence of high outliers for table S5
l<- strsplit(seq, split = "")

# make a matrix of outlier sequences
output <- matrix(unlist(l), ncol = 31, byrow = TRUE)

cbind(pos, strsplit(seq, split = ""))



pos<- c(508994,552251, 587402, 
        701942, 836753,976517,
        1061432,112312,1630026,
        1653648,1828950)
seq<- outlier_lower.seq

# get sequence of low outliers for table S5
l<- strsplit(seq, split = "")

# make a matrix of outlier sequences
output <- matrix(unlist(l), ncol = 31, byrow = TRUE)

