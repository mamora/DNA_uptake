######################################################################
## create a model that predicts DNA uptake in NPlarge frag data #####
######################################################################

# the objective of this script is to run the model to predict uptake


####################################################################################
######             1. load samples and working directory                           #         
library(dplyr)                                                                     #
library(tidyr)                                                                     #
library(foreach)                                                                   #
library(doParallel)                                                                #
                                                                                   #
# uptake ratios                                                                    #
Uptake.ratio.np<- read.csv("./Uptake.ratio.np.corrected.csv") #read uptake file    #
                                                                                   #
# list of USS10                                                                    #
Up.USS.np.10.list<- read.csv("./Uptake.uss.10.list.np.csv") #save file             #
                                                                                   #
# relative number of fragments of input fragment distributions                     #
#as measured by bioanalyzer                                                        #
d.df<- fread("./bioanalyzer_frag_dist.csv")                                        #
####################################################################################


# Settings

settings<- read.csv("./Settings.csv")


###################################################
#                 Settings                        #
###################################################
# set the genomic region to be predicted          #
st = settings$St                                  #
end = settings$end                                #
###################################################


########################################################################################################
#                                            circularize the USS list                                  #
#  Add  the last 10 USS to the beginning of the list and the first 9 USS to the end of the list.       #
#                                                                                                      #
#  Replace the genomic position corresponding to the pos 16 (keypos) of each uss list by               #
#  the corresponding positions if the genome was larger (end) or by the difference of the pos 16       #
#  to the end of the genome (beginning)                                                                #
#                                                                                                      #    
# basically, I will:                                                                                   #
# Add the 10 last uss and the first 9 uss to the uss list                                              #  
                                                                                                       #      
s0<- head(Up.USS.np.10.list$keypos, 9) #first 9 uss                                                    #
s1<- tail(Up.USS.np.10.list$keypos, 10) #last 9 uss                                                    #
                                                                                                       #   
# Replace the genomic position corresponding to the pos 16 (keypos) of each uss list                   #
# by the corresponding positions if the genome was larger (end) or by the difference of the pos 16 to  #
# the end of the genome (beginning)                                                                    #
s2<- (length(Uptake.ratio.np$pos) + s0)                                                                #  
s3<- -(length(Uptake.ratio.np$pos) - s1)                                                               #
circle.uss<- c(s3,Up.USS.np.10.list$keypos, s2)                                                        #
circle.uss.list<- rbind(tail(Up.USS.np.10.list, 10), Up.USS.np.10.list, head(Up.USS.np.10.list, 9) )   #
circle.uss.list$keypos<- circle.uss                                                                    #
########################################################################################################



##################################################################################################
###      From the relative number of fragments, calculate the density distribution               #
# of NP small fragment size and pick 50 - 60 points of the distribution to use in the function.  #

# this section needs to be personalized according to the sample to be analyzed

# pick only NP small fragment distribution
up13<- d.df %>% dplyr::filter(sample == "UP13")
up13.dist<- up13[16:335,]
up13.dist<- up13.dist[up13.dist$bases <= 12000,]

# Fit loess smooth curve to the rel # of fragments
p1 <- ggplot(aes(x = bases, y = rel_mol), data = up13.dist) +
  geom_line() + #span controls the amount of smoothing
  geom_smooth(method = "loess", span = 0.2, n = 10) + #span controls the amount of smoothing
  scale_x_continuous(limits = c(0,16000), breaks = c(0,1000,3000,5000,7000,10000,16000), expand = c(0, 0))+
  ggtitle("Fragment distribution of short input fragments") +
  labs(x = "base pairs", y = "Relative # of fragments") +
  theme_bw() +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))

# save a dataframe with 60 points if fitted loess y-hat data
data.d<- ggplot_build(p1)$data[[2]]

# make sure first point starts at 0
data.d$y[which(data.d$y < 0)] <- 0


# calculate the density of the distribution based on the fitted relative number of fragments
density<- data.d$y/(sum(data.d$y))

sum(density) # should be 1

up13.dist1<- data.table(bases = data.d$x, density = density)

sum(up13.dist1$density)

p1 <- ggplot(aes(x = bases, y = density), data = up13.dist1) +
  geom_point() +
  geom_line() + #span controls the amount of smoothing
  scale_x_continuous(limits = c(0,16000), breaks = seq(0 , 16000, 5000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0,0.4), expand = c(0, 0)) +
  ggtitle("Fragment distribution of short input fragments") +
  labs(x = "base pairs", y = "density") +
  theme() +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))

p1

write.csv(up13.dist1, "./Marcelo_paper_figures_new_order/Figure_6/July_01_2018/up13.dist1.csv")


##########################################################################################################
#                         Use sigmodial equation to calculate probability to be taken up                 #  
                                                                                                         #
low_asy = 0.099  # lower_asymptote                                                                       #
top_asy = 3.751  # top_asymptote (carriyng capacity)                                                     #
grow_rate = 5.5  # growth rate                                                                           #
time_max = 10.631 # time_max (maximum growth point)                                                      #
j<- seq(0, 12.7, by = 0.01) # sequence of USS scores                                                     #
# calculate probability of being taken up based on USS score                                             #
d1<- c()                                                                                                 #
for (i in 1:length(j)){                                                                                  #
  d0<- low_asy + ((top_asy - low_asy)/(1 + exp(-grow_rate * (j[i] - time_max))))                         #
  d1<- c(d1, d0)                                                                                         #
}                                                                                                        #
#adjust to  a range from 0 to 1                                                                          #
d1 <-  d1/top_asy                                                                                        #
                                                                                                         #
p.taken.up<- c()                                                                                         #
for(i in 1:length(circle.uss.list$USS.score)){                                                           #
  pu<- d1[which(abs(j-circle.uss.list$USS.score[i])==min(abs(j- circle.uss.list$USS.score[i])))]         #
  p.taken.up<- c(p.taken.up, pu)                                                                         #
}                                                                                                        #
                                                                                                         #
circle.uss.list$p.taken.up<- p.taken.up                                                                  #
##########################################################################################################




#####################################################################################################################
#       incorporate the model function to a loop that calculates contribution to uptake                             #
#                               of all fragment sizes for all genomic positions                                     #  
source("./model_pred_up.R")                                                                                         #
                                                                                                                    #
cl<-makeCluster(4)                                                                                                  #
registerDoParallel(cl)                                                                                              #
                                                                                                                    #
run.d<-  function(a){                                                                                               #
  h<- foreach(ds = 1:length(up13.dist1$bases), .combine='c') %do% {                                                 #
    prob.model_f(a = a, b = circle.uss.list$keypos, w = up13.dist1$bases[ds], fr = up13.dist1$density[ds], sd = sd) #
  }                                                                                                                 #
  return(h)                                                                                                         #
}                                                                                                                   #
                                                                                                                    #
d<- foreach(a = st:end, .combine='cbind',.packages= 'doParallel') %dopar% {                                         #
  run.d(a)                                                                                                          #
}                                                                                                                   #
write.csv(d, file = "test_20kb_np.csv")                                                                             #
#####################################################################################################################



