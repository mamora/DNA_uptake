###################################################################################################
#      analyze the effect of USS density and of USS score on uptake of large fragments data       #
###################################################################################################

# load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)
library(zoo)

# load datasets

folder.name<- "./Marcelo_paper_figures_new_order/large_frag_analysis/"
# load scores 
np.USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")
# load uptake ratios
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv") #read uptake file 
# list of USS10
Up.USS.np.10.list <- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) 

#####################################################################
#              define a window to do the analysis                   #
#####################################################################

#  To improve speed of functions and loops
window<- 3001  # window size, change this as see fit

t.n<- c(rep(NA,length(np.USS.scores$V1)))  # empty vector of NP genome size
g.n<- seq(1,length(t.n)) #vector for loop

#####################################################################
#    define functions to find USS and their score in a fragment     #
#####################################################################

#function to find matches between USS and a defined section
find<- function(x){
  length(intersect(x, Up.USS.np.10.list$USS.pos))
}


# function to get mean score of uss in a window
find.m.score.np<- function(x){
  f<- intersect(x, Up.USS.np.10.list$USS.pos)
  mean(np.USS.scores$max[f])
}

# function to get mean uptake ratio of a USS in a window
find.m.ratio.np<- function(x){
  mean(Uptake.ratio.np$ratio_long[x])
}


#####################################################################
#                     Run the  analysis                             #
#####################################################################

# get NP positions circularized
seq1<- c(np.USS.scores$V1[(length(np.USS.scores$V1) - ((window/2)-1)):length(np.USS.scores$V1)], np.USS.scores$V1, np.USS.scores$V1[1:((window/2)-1)])  # join sequence with first 100 bases to circularize

# get number of USS per 3kb in sliding window
den_uss_3kb.np<- rollapply(seq1, by = 1, width = window, FUN = find, align = "center", partial = FALSE)  # run function centered with sliding window for GG
# get the score of those USS
den_uss_3kb.np.score<- rollapply(seq1, by = 1, width = window, FUN = find.m.score.np, align = "center", partial = FALSE)  # run function centered with sliding window for GG
# get the uptake ratio of those USS
den_uss_3kb.np.ratio<- rollapply(seq1, by = 1, width = window, FUN = find.m.ratio.np, align = "center", partial = FALSE)  # run function centered with sliding window for GG

# put everything in a dataframe
density_3kb<- data.table(density = den_uss_3kb.np, mean.score = den_uss_3kb.np.score,
                         mean.ratio = den_uss_3kb.np.ratio)

write.csv(density_3kb, "./Marcelo_paper_figures_new_order/large_frag_analysis/den_3kb.csv⁩" )

density_3kb$score.f<- rep("NA", length(density_3kb$density))

# subset by mean score
density_3kb$score.f[which(density_3kb$mean.score >= 10 & density_3kb$mean.score < 10.5)] <- ">= 10 < 10.5"
density_3kb$score.f[which(density_3kb$mean.score >= 10.5 & density_3kb$mean.score < 11)] <- ">= 10.5 < 11"
density_3kb$score.f[which(density_3kb$mean.score >= 11 & density_3kb$mean.score < 11.5)] <- ">= 11 < 11.5"
density_3kb$score.f[which(density_3kb$mean.score >= 11.5)] <- ">= 11.5"

str(density_3kb)

# factorize variables
density_3kb$score.f<- as.factor(density_3kb$score.f)
density_3kb$density<- as.factor(density_3kb$density)

# plot density over score of 3kb fragments
p<- density_3kb  %>%    
  ggplot(aes(x = density, y = mean.ratio, fill =  score.f)) +
  geom_boxplot() +
  scale_x_discrete() +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))+
  labs(y = "mean uptake ratio", x = "number of USS") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=24),
        axis.title = element_text(size = 24, face = "bold")) 

p

file_name = paste(folder.name,"den_3kb","NP", ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, units = "px")
print(p)
dev.off()