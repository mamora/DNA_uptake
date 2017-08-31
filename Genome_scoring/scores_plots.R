
####################     Map of USS scores vs locations      #########################
######################################################################################

# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/" # Where these files are located

# Set working directory
setwd(whereami)

list.files(whereami) #see files in directory

#############################load packages##############################
library(ggplot2)
library(dplyr)
library(zoo)
library(data.table)
#########################################################################

np.uptake.scores<- read.csv("./datasets/np.uptake.scores.csv") #read scores file 

gg.uptake.scores<- read.csv("./datasets/PittGG.uptake.scores.csv") #read scores file 

scores.with.density.np<- read.csv("./datasets/scores.with.density.np.csv")

str(np.uptake.scores)


np.uptake.scores$X <- as.numeric(np.uptake.scores$X)


u<- np.uptake.scores %>%  dplyr::slice(1:10000) %>% 
  ggplot(aes(x = X, y = c)) +
  geom_point(aes(colour = c),shape = 20, size = 1)+
  scale_colour_gradient(low = "white", high = "black", name = "uptake score") +
  scale_x_continuous(breaks = seq(0 , 10000, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 13), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake scores") +
  ggtitle(" Uptake scores vs locations in the NP genome") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/ratiovslocations.np","10kb", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(u)
dev.off()


u<- np.uptake.scores  %>% 
  ggplot(aes(x = X, y = c)) +
  geom_point(aes(colour = c),shape = 20, size = 1)+
  scale_colour_gradient(low = "white", high = "black", name = "uptake score") +
  scale_x_continuous(breaks = seq(0 , 1914386, 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(8, 13), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "uptake scores") +
  ggtitle(" Uptake scores vs locations in the NP genome") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/ratiovslocations.np","whole_genome","only_high", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(u)
dev.off()

##  mmmm I do not see this figure very useful


####################################################################################


Up.USS.np.10.list<- read.csv("./datasets/Up.USS.np.10.list.csv")

Up.USS.PittGG.10.list<- read.csv("./datasets/Up.USS.PittGG.10.list.csv")

length(which(np.uptake.scores$X[1:20000] %in% Up.USS.np.10.list$USS.pos)) #  this script finds all positions that have a uss. But if 1 position have 2 uss's it only accounts for one. 
 
length(which(Up.USS.np.10.list$USS.pos %in% np.uptake.scores$X[1:20000])) #  this script finds all uss that are in a range of positions. if 1 position have 2 uss's it accounts both of them.

window<- 2000  # window size

t.n<- c(rep(NA,length(np.uptake.scores$X)))  # empty vector of NP genome size
t.g<- c(rep(NA,length(gg.uptake.scores$X)))  # empty vector of PittGG genome size
 
data.np<- c(np.uptake.scores$X, np.uptake.scores$X[1:window])

data.gg<- c(gg.uptake.scores$X, gg.uptake.scores$X[1:window.size])

g.n<- seq(1,length(t.n)) #vector for loop

g.g<- seq(1,length(t.g)) #vector for loop


i = 1

system.time(for (i in g){
  s<- 0 + i
  e<- window + i  
  t[i]<-  length(which(data.np[s:e] %in% Up.USS.np.10.list$USS.pos))
})

find(data.np[s:e])

#function to find matches between USS and a defined section of numbers
find<- function(x){
  length(intersect(x, Up.USS.np.10.list$USS.pos))
}

# circularize number of positions
seq1<- c(np.uptake.scores$X[(length(np.uptake.scores$X) - ((window/2)-1)):length(np.uptake.scores$X)], np.uptake.scores$X, np.uptake.scores$X[1:((window/2)-1)])  # join sequence with first 100 bases to circularize


cg1<- rollapply(seq1, by = 1, width = window, FUN = find, align = "center", partial = FALSE)



np.uptake.scores$t.2kb<- cg1


write.csv(np.uptake.scores, "./datasets/scores.with.density.np.csv")







t1<- c(rep(NA,length(gg.uptake.scores$X)))

data.s1<- c(gg.uptake.scores$X, gg.uptake.scores$X[1:20000])

length(data)

g1<- seq(1,length(t1))



system.time(for (i in g1){
  s<- 0 + i
  e<- 20000 + i  
  t[i]<-  length(which(data.s1[s:e] %in% Up.USS.PittGG.10.list$USS.pos))
})




gg.uptake.scores$density<- t1


write.csv(np.uptake.scores, "./datasets/scores.with.density.np.csv")






str(scores.with.density.np)


p<- scores.with.density.np %>% 
  ggplot(aes(x = X, y = t)) +
  geom_line(size = 1)+
  scale_x_continuous(breaks = seq(0, length(scores.with.density.np$X), 200000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0))+
  labs(x = "NP genome positions", y = "USS density") +
  ggtitle("NP USS density over a 20kb window") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold")) 
file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/USS_density","20kb_wind","NP", "whole_genome", "tiff", sep=".")
tiff(file_name, width = 1000, height = 700, units = "px")
print(p)
dev.off() 


########################################      test statistically distributions of USS scores  ################################

# I will use a Kolmogorov-Smirnov test which is a non-parametric test (does not assume normal distribution) 
# that compares the media, deviation and  the shape of the distributions to be tested.


ks.test(np.uptake.scores$c, gg.uptake.scores$c)

ks.test(np.uptake.scores$w, gg.uptake.scores$w)







