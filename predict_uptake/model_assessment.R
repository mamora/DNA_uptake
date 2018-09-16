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
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv") #read uptake file 

# list of USS10
Up.USS.np.10.list<- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) #save file

up15.dist1<- fread("./Marcelo_paper_figures_new_order/model/Sept_10_2018/up15.dist_adj_4.w.csv")


pred_100kb_w4<- fread("./Marcelo_paper_figures_new_order/model/Sept_10_2018/new_model_100kb_np_dist_w4.csv")

pred_100kb_w4$V1<- NULL

# sum contribution of all fragment alignments
total.p1<- apply(pred_100kb_w4, MARGIN = 1, sum)


sum(pred_100kb_w4[1,])

sum(pred_100kb_w4[100000,])


# for NP
st = 1   #
end = 100000 


mean(total.p1) 


mean(total.p1[1:100000])



mean(Uptake.ratio.np$ratio_short[1:100000])


#######################################################
# 5. plot real data + predicted data in the same plot #
#######################################################

np.USS.scores<- fread("./datasets/final_datasets/USS_scores/uptake_model/np.USS.scores.corrected.csv")

###################################
#   Normalize predicted uptake    #
###################################
norm<- function (data = model$pred.p.0.001){
  h_pe<- mean(data, na.rm = FALSE)
  s_pe<- (data * 1)/h_pe
  return(s_pe)
}

data = total.p1

mean(data, na.rm = FALSE)


s1<- norm(data = total.p1)


ob<- norm(data = Uptake.ratio.np$ratio_short[st:end])


mean(s1)
mean(ob)


# function to add arrows in places where there are USS

make.arrows.plot<- function(p0, start, endf){
  data1<- np.USS.scores[(start):(endf),]
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
    y = 0.2, #seq(0.2, (0.2 *length(st)),0.2),
    yend = 0.2 #seq(0.2, (0.2 *length(st)),0.2)
  )
  
  p0<- p0 + geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), colour = "red")
  return(p0)
}


real.data<- data.table(pos = Uptake.ratio.np$pos[1:20000], observed_uptake = ob[1:20000], expected_uptake = s1[1:20000])

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


file_name = paste("./Marcelo_paper_figures_new_order/model/Sept_10_2018/new_dist_100kb_4","NP", "tiff", sep=".")
tiff(file_name, width = 1200, height = 800, units = "px")
print(p1)
dev.off()



distributions<- rbind(up15.dist_reads, up15.dist_bioanalyzer)
dist.names<- c(rep("up15.dist_reads", 58), rep("up15.dist_bioanalyzer", 60))
distributions$names<- dist.names

distributions$bases<- as.numeric(distributions$bases)
 
p<-  ggplot(aes(x = bases, y = density), data = distributions) +
  geom_point() +
  geom_line() +
  facet_grid(names~.) +
  labs(x = "base pairs", y = "density") +
  ggtitle("distributions used to model predicted uptake") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=22),
        axis.title = element_text(size = 22, face = "bold")) 


file_name = paste("./Marcelo_paper_figures_new_order/model/July_18_2018/dist_hist","NP", "tiff", sep=".")
tiff(file_name, width = 1200, height = 800, units = "px")
print(p)
dev.off()



