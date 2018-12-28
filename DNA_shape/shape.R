###########################################
##        DNA shape analysis             ##
###########################################


# load packages
library(dplyr)
library(plyr)
library(tidyr)
library(data.table)
library(DNAshapeR)
library(ggplot2)
library(ggpubr)


# as if the whole file were copy-pasted to the R command-line


########################
#   load datasets      #
########################
folder.name <- "./Marcelo_paper_figures_new_order/DNA_shape/"  
Up.USS.np.9.5.list<- fread(here::here("datasets/final_datasets","Up.USS.np.9.5.list.csv")) 
# Name the sites file, and genome reference file
fasta.np    <- "datasets/sequences/np.pilon.fasta" # Reference fasta to the NP genome
# Read in genome
genome.np <- read.genome(fasta.np) # read fasta with seqinr, then massage

########################
#   load functions     #
########################
# Load the various functions used below from the "pssmFunctions.R" file,
# as if the whole file were copy-pasted to the R command-line
source("./helper_functions/pssmFunctions1.R")


# function to extract the sequence of USS 
make.seq.list<- function(low.uss, save){
  for (i in 1:length(low.uss$USS.pos)){
    if (low.uss$strand[i] == "w") { 
      site = low.uss$USS.pos[i]
      w<- as.character(genome.np[[1]][(site - 2):(site + 32)])
      write.fasta(sequences = w, names = low.uss$USS.pos[i], open = "a", file.out = save)
      low.uss.10<- c(low.uss.10,paste(w, collapse = ""))
    } else {
      site = low.uss$USS.pos[i]
      g<- genome.np[[1]][(site + 32):(site - 2)]
      c<-  as.character(revalue(g, c("a"="t", "t"="a","g"="c","c"="g","n" = "n")))
      write.fasta(sequences = c, names = low.uss$USS.pos[i], open = "a", file.out = save)
      low.uss.10<- c(low.uss.10,paste(c, collapse = ""))
    }  
  }  
  return(low.uss.10)
}

###############################################################################
#   subset isolated USS by uptake ratio level  and extract their sequence     #
###############################################################################

# strong USS
high<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$next.uss >= 800 & Up.USS.np.9.5.list$keyup_small >= 3),]
# medium USS
mid<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$next.uss >= 800 & Up.USS.np.9.5.list$keyup_small >= 0.5 & Up.USS.np.9.5.list$keyup_small < 3),]
# weak USS
low<- Up.USS.np.9.5.list[which(Up.USS.np.9.5.list$next.uss >= 800 & Up.USS.np.9.5.list$keyup_small < 0.5),]

# extract sequences
low.uss.10<- NULL
save<- paste(folder.name,"all_uss",".fasta", sep = "")
uss.h.seq<- make.seq.list(Up.USS.np.9.5.list, save)

low.uss.10<- NULL
save<- paste(folder.name,"high",".fasta", sep = "")
uss.h.seq<- make.seq.list(high, save)

low.uss.10<- NULL
save<- paste(folder.name,"mid",".fasta", sep = "")
uss.h.seq<- make.seq.list(mid, save)

low.uss.10<- NULL
save<- paste(folder.name,"low",".fasta", sep = "")
uss.h.seq<- make.seq.list(low, save)

################################################
#  analyze USS DNA shape structural features   #
################################################

#  analyze DNA shapes
high.s<- getShape("./Marcelo_paper_figures_new_order/DNA_shape/high.fasta")
mid.s<- getShape("./Marcelo_paper_figures_new_order/DNA_shape/mid.fasta")
low.s<- getShape("./Marcelo_paper_figures_new_order/DNA_shape/low.fasta")



# plot minimum groove width 
plotShape(high.s$MGW)
plotShape(mid.s$MGW)
plotShape(low.s$MGW)

# plot helix twist
plotShape(high.s$HelT)
plotShape(mid.s$HelT)
plotShape(low.s$HelT)

# plot Propeller twist
plotShape(high.s$ProT)
plotShape(mid.s$ProT)
plotShape(low.s$ProT)

# plot Base roll
plotShape(high.s$Roll)
plotShape(mid.s$Roll)
plotShape(low.s$Roll)

#############################################################
#   analyze minimum groove width of USS with different ratios
#############################################################

# make a long format dataframe of MGW values  
high.s1.mgw <- stack(as.data.frame(t(high.s$MGW[,3:33])))
mid.s1.mgw <- stack(as.data.frame(t(mid.s$MGW[,3:33])))
low.s1.mgw <- stack(as.data.frame(t(low.s$MGW[,3:33])))

# add which USS have high, medium or low ratios
high.s1.mgw$ratio <- rep("high", length(high.s1.mgw$values))
mid.s1.mgw$ratio <- rep("mid", length(mid.s1.mgw$values))
low.s1.mgw$ratio <- rep("low", length(low.s1.mgw$values))

# add USS sequences positions to dataframe
high.s1.mgw$seq.pos<- rep(1:31, length(unique(high.s1.mgw$ind)))
mid.s1.mgw$seq.pos<- rep(1:31, length(unique(mid.s1.mgw$ind)))
low.s1.mgw$seq.pos<- rep(1:31, length(unique(low.s1.mgw$ind)))

MGW.m<- rbind(high.s1.mgw,mid.s1.mgw,low.s1.mgw)
MGW.m$ind<-  as.factor(MGW.m$ind) 

# plot MGW DNA shape over USS positions, subset by uptake ratio
p1<- MGW.m %>% 
  ggplot(aes(x = seq.pos, y =values)) +
  geom_line(aes(group = ind), size = 1) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 8), expand = c(0, 0)) +
  facet_grid(ratio ~.) +
  labs(y = "minor groove width (MWM) angstrom", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 30),
        axis.title = element_text(size = 30, face = "bold")) 

p1 

file_name = paste("./outliers/test",i,"NP", "tiff", sep =".")
tiff(file_name, width = 1200, height = 1000, units = "px")
print(pt)
dev.off()


#############################################################
#   analyze helix twist width of USS with different ratios  #
#############################################################

# make a long format dataframe of HelT values  
high.s1.HelT <- stack(as.data.frame(t(high.s$HelT[,3:33])))
mid.s1.HelT <- stack(as.data.frame(t(mid.s$HelT[,3:33])))
low.s1.HelT <- stack(as.data.frame(t(low.s$HelT[,3:33])))

# add which USS have high, medium or low ratios
high.s1.HelT$ratio <- rep("high", length(high.s1.HelT$values))
mid.s1.HelT$ratio <- rep("mid", length(mid.s1.HelT$values))
low.s1.HelT$ratio <- rep("low", length(low.s1.HelT$values))

# add USS sequences positions to dataframe
high.s1.HelT$seq.pos<- rep(1:31, length(unique(high.s1.HelT$ind)))
mid.s1.HelT$seq.pos<- rep(1:31, length(unique(mid.s1.HelT$ind)))
low.s1.HelT$seq.pos<- rep(1:31, length(unique(low.s1.HelT$ind)))

HelT.m<- rbind(high.s1.HelT,mid.s1.HelT,low.s1.HelT)

HelT.m$ind<-  as.factor(HelT.m$ind) 


# plot HelT DNA shape over USS positions, subset by uptake ratio
p1<- HelT.m %>% 
  ggplot(aes(x = seq.pos, y =values)) +
  geom_line(aes(group = ind), size = 1) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(30, 40), expand = c(0, 0)) +
  facet_grid(ratio ~.)  +
  labs(y = "helix twist (degree)", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 10),
        axis.title = element_text(size = 10, face = "bold")) 

p1 

#################################################################
#   analyze Propeller twist width of USS with different ratios  #
#################################################################

# make a long format dataframe of ProT values  
high.s1.ProT <- stack(as.data.frame(t(high.s$ProT[,3:33])))
mid.s1.ProT <- stack(as.data.frame(t(mid.s$ProT[,3:33])))
low.s1.ProT <- stack(as.data.frame(t(low.s$ProT[,3:33])))

# add which USS have high, medium or low ratios
high.s1.ProT$ratio <- rep("high", length(high.s1.ProT$values))
mid.s1.ProT$ratio <- rep("mid", length(mid.s1.ProT$values))
low.s1.ProT$ratio <- rep("low", length(low.s1.ProT$values))

# add USS sequences positions to dataframe
high.s1.ProT$seq.pos<- rep(1:31, length(unique(high.s1.ProT$ind)))
mid.s1.ProT$seq.pos<- rep(1:31, length(unique(mid.s1.ProT$ind)))
low.s1.ProT$seq.pos<- rep(1:31, length(unique(low.s1.ProT$ind)))

ProT.m<- rbind(high.s1.ProT,mid.s1.ProT,low.s1.ProT)

ProT.m$ind<-  as.factor(ProT.m$ind) 

# plot ProT DNA shape over USS positions, subset by uptake ratio
p1<- ProT.m %>% 
  ggplot(aes(x = seq.pos, y =values)) +
  geom_line(aes(group = ind), size = 1) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-20, 5), expand = c(0, 0)) +
  facet_grid(ratio ~.)  +
  labs(y = "propeller twist (degree)", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 10),
        axis.title = element_text(size = 10, face = "bold")) 

p1


#############################################################
#   analyze Base Roll width of USS with different ratios    #
#############################################################

# make a long format dataframe of Roll values  
high.s1.Roll <- stack(as.data.frame(t(high.s$Roll[,3:33])))
mid.s1.Roll <- stack(as.data.frame(t(mid.s$Roll[,3:33])))
low.s1.Roll <- stack(as.data.frame(t(low.s$Roll[,3:33])))

# add which USS have high, medium or low ratios
high.s1.Roll$ratio <- rep("high", length(high.s1.Roll$values))
mid.s1.Roll$ratio <- rep("mid", length(mid.s1.Roll$values))
low.s1.Roll$ratio <- rep("low", length(low.s1.Roll$values))

# add USS sequences positions to dataframe
high.s1.Roll$seq.pos <- rep(1:31, length(unique(high.s1.Roll$ind)))
mid.s1.Roll$seq.pos <- rep(1:31, length(unique(mid.s1.Roll$ind)))
low.s1.Roll$seq.pos <- rep(1:31, length(unique(low.s1.Roll$ind)))

Roll.m <- rbind(high.s1.Roll,mid.s1.Roll,low.s1.Roll)

Roll.m$ind <-  as.factor(Roll.m$ind) 

# plot base roll DNA shape over USS positions, subset by uptake ratio
p1 <- Roll.m %>% 
  ggplot(aes(x = seq.pos, y =values)) +
  geom_line(aes(group = ind), size = 1) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-15, 10), expand = c(0, 0)) +
  facet_grid(ratio ~.)  +
  labs(y = "Roll (degree)", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 10),
        axis.title = element_text(size = 10, face = "bold")) 

p1

##############################################################################
#       calculate mean and standard deviation for each structural feature    #
##############################################################################

#####################
#       MGW         #
#####################
mean1 <- MGW.m %>% filter(ratio == "high") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean2 <- MGW.m %>% filter(ratio == "mid") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean3 <- MGW.m %>% filter(ratio == "low") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
sd1 <- MGW.m %>% filter(ratio == "high") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd2 <- MGW.m %>% filter(ratio == "mid") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd3 <- MGW.m %>% filter(ratio == "low") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))

pos <- c(mean1$seq.pos, mean2$seq.pos, mean3$seq.pos)
mean<- c(mean1$values, mean2$values,mean3$values)
sd<- c(sd1$values, sd2$values, sd3$values)

data.h <- data.table(pos, mean, sd)
data.h$ratio <- c(rep("high", 31), rep("mid", 31), rep("low", 31)) 

p1<- data.h %>% 
  ggplot(aes(x = pos, y =mean)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 7), expand = c(0, 0)) +
  facet_grid(ratio ~.) +
  labs(y = "minor groove width (MWM) angstrom", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 10),
        axis.title = element_text(size = 10, face = "bold")) 

p1

p1<- data.h %>% 
  ggplot(aes(x = pos, y =mean, colour = ratio)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 7), expand = c(0, 0)) +
  labs(y = "minor groove width (MWM) angstrom", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 20),
        axis.title = element_text(size = 20, face = "bold")) 

p1

#####################
#       ProT        #
#####################

mean1<- ProT.m %>% filter(ratio == "high") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean2<- ProT.m %>% filter(ratio == "mid") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean3<- ProT.m %>% filter(ratio == "low") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
sd1<- ProT.m %>% filter(ratio == "high") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd2<- ProT.m %>% filter(ratio == "mid") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd3<- ProT.m %>% filter(ratio == "low") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))

pos <- c(mean1$seq.pos, mean2$seq.pos, mean3$seq.pos)
mean<- c(mean1$values, mean2$values,mean3$values)
sd<- c(sd1$values, sd2$values, sd3$values)

data.h <- data.table(pos, mean, sd)
data.h$ratio <- c(rep("high", 31), rep("mid", 31), rep("low", 31)) 

p2<- data.h %>% 
  ggplot(aes(x = pos, y =mean, colour = ratio)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-20, 0), expand = c(0, 0)) +
  labs(y = "propeller twist (degree)", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 20),
        axis.title = element_text(size = 20, face = "bold")) 

p2

#####################
#       HelT        #
#####################

mean1<- HelT.m %>% filter(ratio == "high") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean2<- HelT.m %>% filter(ratio == "mid") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean3<- HelT.m %>% filter(ratio == "low") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
sd1<- HelT.m %>% filter(ratio == "high") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd2<- HelT.m %>% filter(ratio == "mid") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd3<- HelT.m %>% filter(ratio == "low") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))

pos <- c(mean1$seq.pos, mean2$seq.pos, mean3$seq.pos)
mean<- c(mean1$values, mean2$values,mean3$values)
sd<- c(sd1$values, sd2$values, sd3$values)

data.h <- data.table(pos, mean, sd)
data.h$ratio <- c(rep("high", 31), rep("mid", 31), rep("low", 31)) 

p3<- data.h %>% 
  ggplot(aes(x = pos, y =mean, colour = ratio)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(30, 40), expand = c(0, 0)) +
  labs(y = "helix twist (degree)", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 20),
        axis.title = element_text(size = 20, face = "bold")) 

p3

#####################
#       Roll        #
#####################

mean1<- Roll.m %>% filter(ratio == "high") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean2<- Roll.m %>% filter(ratio == "mid") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
mean3<- Roll.m %>% filter(ratio == "low") %>%  group_by(seq.pos) %>% summarise_at(vars(values), funs(mean(., na.rm=TRUE)))
sd1<- Roll.m %>% filter(ratio == "high") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd2<- Roll.m %>% filter(ratio == "mid") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))
sd3<- Roll.m %>% filter(ratio == "low") %>% group_by(seq.pos) %>% summarise_at(vars(values), funs(sd(., na.rm=TRUE)))

pos <- c(mean1$seq.pos, mean2$seq.pos, mean3$seq.pos)
mean<- c(mean1$values, mean2$values,mean3$values)
sd<- c(sd1$values, sd2$values, sd3$values)

data.h <- data.table(pos, mean, sd)
data.h$ratio <- c(rep("high", 31), rep("mid", 31), rep("low", 31)) 

p4<- data.h %>% 
  ggplot(aes(x = pos, y =mean, colour = ratio)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  scale_x_continuous(breaks = seq(1,31, by = 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  labs(y = "Roll (degree)", x = "bases of the USS motif") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "right",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= 20),
        axis.title = element_text(size = 20, face = "bold")) 

p4

p5<- ggarrange(p1, p2, p3, p4, 
          ncol = 2, nrow = 2)

file_name = paste(folder.name,"NP_shapes", ".tiff", sep ="")
tiff(file_name, width = 1500, height = 1000, units = "px")
print(p5)
dev.off()


#########################################################################
###                   stats run U mann whitney                         ##
#########################################################################

str(MGW.m)
str(HelT.m)
str(Roll.m)
str(ProT.m)


# This function calculates Mann Whitney U test to each structural feature USS list
sig.pos.f<- function(model = MGW.m){
  o<- c()
  model$ratio<- as.factor(model$ratio)
  model<- model %>% filter(ratio == c("high", "low"))
  model$ratio<- as.numeric(model$ratio)
  
  for(i in 1:length(unique(model$seq.pos))){
    test<- model %>% filter(seq.pos == i)
    y = test$values
    x = test$ratio
    st<- kruskal.test(y ~ x) 
    o<- c(o, st$p.value)
  }
  sig.pos<- which(o < 0.01) 
  return(sig.pos)
}

sig.pos.f(model = MGW.m)
sig.pos.f(model = HelT.m)
sig.pos.f(model = Roll.m)
sig.pos.f(model = ProT.m)



