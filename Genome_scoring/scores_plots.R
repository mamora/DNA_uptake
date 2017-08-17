
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
library(data.table)
#########################################################################

np.uptake.scores<- read.csv("./datasets/np.uptake.scores.csv") #read scores file 

gg.uptake.scores<- read.csv("./datasets/PittGG.uptake.scores.csv") #read scores file 

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

t<- c(rep(NA,length(np.uptake.scores$X)))

number.uss<- function(data = np.uptake.scores$X, uss.pos = Up.USS.np.10.list$USS.pos){
num.uss<- length(which(uss.pos %in% data)) #  this script finds all uss that are in a range of positions. if 1 position have 2 uss's it accounts both of them.
}

 
data.s<- c(np.uptake.scores$X, np.uptake.scores$X[1:20000])

length(data)

g<- seq(1,length(t))

g1<- seq(1:1000)

system.time(for (i in 1:1000){
  s<- 0 + i
  e<- 20000 + i  
  t[i]<-  length(which(data.s[s:e] %in% Up.USS.np.10.list$USS.pos))
})


system.time(for (i in g1){
  s<- 0 + i
  e<- 20000 + i  
  t[i]<-  length(which(data.s[s:e] %in% Up.USS.np.10.list$USS.pos))
})

system.time(for (i in g){
  s<- 0 + i
  e<- 20000 + i  
  t[i]<-  length(which(data.s[s:e] %in% Up.USS.np.10.list$USS.pos))
})




str(t)

f<-  number.uss(data = np.uptake.scores$X[1:20002], uss.pos = Up.USS.np.10.list$USS.pos)


library(microbenchmark)





# works but too slow
d<- rollapply(np.uptake.scores$X[1:20001], width = 20000, by = 1, align = "left", FUN = number.uss, uss.pos = Up.USS.np.10.list$USS.pos )


