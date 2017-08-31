# Name the path of the working directory; 
### Set the path from your computer ###
whereami <- "C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/" # Where these files are located

# Set working directory
setwd(whereami)

library(ggplot2)
library(hexbin)
library(RColorBrewer)

Uptake.ratio.np<- read.csv("./datasets/Uptake.ratio.np.csv") #read scores file 

Uptake.ratio.gg<- read.csv("./datasets/Uptake.ratio.gg.csv") #read scores file 

Uptake.uss.10.list.np<- read.csv("./datasets/Up.USS.np.10.list.csv")



#plot the uptake ratio of short fragments vs distance to USS

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

test<- Uptake.ratio.np$pos + Uptake.ratio.np$close.USS.np


test1<- match(test, Uptake.uss.10.list.np$keypos)

g<- which(is.na(test1 == TRUE)) 


PittGG.uptake.ratio$sim_dist_uss<- PittGG.uptake.ratio$close.USS


#make distance to uss 
PittGG.uptake.ratio$sim_dist_uss[g]<- PittGG.uptake.ratio$close.USS[g] - (PittGG.uptake.ratio$close.USS[g] * 2)



sum(PittGG.uptake.ratio$sim_dist_uss[g] + PittGG.uptake.ratio$close.USS[g]) #did worked should sum to 0


p<- hexbinplot(ratio_short ~ close.USS.np, data=Uptake.ratio.np, 
           colramp=rf, xbins=100, 
           trans=log, inv=exp, aspect=1,
           xlab="Distance from USS central position 16",
           ylab="Uptake Ratio",
           main = "Uptake ratio of NP from small fragments vs distance to USS"
)

file_name = paste("C:/Users/marcelo/Dropbox/uptake/Uptake_summer2017/Figures/ratio_vs_distance_uss","NP", "tiff", sep=".")
tiff(file_name, width = 1400, height = 800, units = "px")
print(p)
dev.off()



hexbinplot(re.ratio_long ~ sim_dist_uss, data=PittGG.uptake.ratio, 
           colramp=rf, xbins=100, 
           trans=log, inv=exp, aspect=1,
           xlab="Distance from USS central position 17",
           ylab="Uptake Ratio",
           main = "Uptake ratio of PittGG from large fragments vs distance to USS"
)

