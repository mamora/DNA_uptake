######################################################################
## create a model that predicts DNA uptake in NP small frag data #####
######################################################################


###########################################################
## Warning: Please open Uptake_summer2017.Rproj first #####
###########################################################

##############################################################
######   1. load samples and working directory        ########         
##############################################################


library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

# uptake ratios                                                                
folder.name <- "./Marcelo_paper_figures_new_order/model/Sept_24_2018/"  
#read uptake file
Uptake.ratio.np <- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")      
# load USS scores
USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")
# list of USS10
Up.USS.np.10.list <- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) 
# fragment sizes
up15.dist1 <- fread(paste(folder.name,"up15.dist_raw_den_class.csv", sep = ""))
# table of predicted uptake calculations
pred_100kb_raw<- fread(paste(folder.name,"simple_model_100kb_np_raw_den_class_2.csv", sep = ""))

# remove extra columns
pred_100kb_raw$V1<- NULL

# sum contribution of all fragment alignments
total.p1<- apply(pred_100kb_raw, MARGIN = 1, sum)

# define genomic segment 
st = 1
end = 100000 


###################################
#   Normalize predicted uptake    #
###################################
# function for normalizing predicted uptake
norm<- function (data = model$pred.p.0.001){
  h_pe<- mean(data, na.rm = FALSE)
  s_pe<- (data * 1)/h_pe
  return(s_pe)
}

# normalize predicted uptake to a mean of 1
s1<- norm(data = total.p1)
# normalize predicted uptake to a mean of 1
ob<- norm(data = Uptake.ratio.np$ratio_short[st:end])

# check if it was normalize correctly. Mean should be 1
mean(s1)
mean(ob)

#######################################################
#  plot real data + predicted data in the same plot #
#######################################################

# function to add arrows in places where there are USS in a plot
make.arrows.plot<- function(p0, start, endf){
  data1<- USS.scores[(start):(endf),]
  f<- which(data1$max >= 10)
  t2 <- data1[f,]
  st<- c()
  end<- c()
  for (j in 1:length(t2$w)){
    if (t2$strand[j] == "w"){
      st[j]<- t2$V1[j]
      end[j]<- t2$V1[j] + 30
    }else{
      st[j]<- t2$V1[j] + 29
      end[j]<- t2$V1[j] - 1
    }
  }
  segment_data = data.frame(
    x = st,
    xend = end, 
    y = 0.2, # place where I want arrows to be in the y axis
    yend = 0.2 
  )
  p0<- p0 + geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),
                         arrow = arrow(length = unit(0.2, "cm"), type = "closed"), colour = "red")
  return(p0)
}

# put observed and predicted uptake in a dataframe
real.data<- data.table(pos = Uptake.ratio.np$pos[1:20000], observed_uptake = ob[1:20000],
                       expected_uptake = s1[1:20000])

# reshape dataframe to long format and then plot uptake
p<- real.data %>% tidyr::gather("sample", "uptake", 2:3) %>%    
  ggplot() +
  geom_point(aes(x = pos, y = uptake, colour = sample ), shape = 20, size = 1) +
  scale_x_continuous(breaks = seq(st , end, 2000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0))+
  labs(x = "NP genomic positions", y = "predicted and observed uptake ratios") +
  ggtitle("predicted small fragment NP uptake ratios with baselines ") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=22),
        axis.title = element_text(size = 22, face = "bold")) 

# add arrows
p1<- make.arrows.plot(p0 = p, start =1, endf = 20000)


file_name = paste(folder.name,"new_dist_20kb_den_class_raw_2","NP", ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, units = "px")
print(p1)
dev.off()




#################################################################
#  plot all positions observed vs predicted uptake Figure 5     #
#################################################################

real.data<- data.table(pos = Uptake.ratio.np$pos[1:100000], observed_uptake = ob[1:100000], 
                       expected_uptake = s1)

p<- real.data  %>%    
  ggplot(aes(x = observed_uptake, y = expected_uptake)) +
  geom_hex(bins = 300) +
  scale_fill_gradientn(trans="log10", colours = rainbow(7))+
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  labs(y = "predicted uptake ratios", x = "observed uptake ratios") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=22),
        axis.title = element_text(size = 22, face = "bold")) 

file_name = paste(folder.name,"obs_vs_pred_100kb_small_den_class","NP", ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, units = "px")
print(p)
dev.off()

#################################################################
#  plot residuals of observed - expected uptake                 #
#################################################################

residuals<- real.data$observed_uptake - real.data$expected_uptake
real.data$residuals<- residuals

p<-  ggplot() +
  geom_line(aes(x = pos, y = residuals), data = real.data) +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "NP genomic positions", y = "residuals") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=22),
        axis.title = element_text(size = 22, face = "bold")) 

p

file_name = paste(folder.name,"residuals_100kb_raw_small","NP", ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, units = "px")
print(p1)
dev.off()


