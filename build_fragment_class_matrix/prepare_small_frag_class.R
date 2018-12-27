##################################################################################################
#     this script uses fragment distribution of large fragment size data measured by bioanalyzer #
#                          then assign fragment sizes into 200bp classes                         #  
##################################################################################################


##########################
#   load the packages    #
##########################
library(data.table)      #
library(ggplot2)         #
library(dplyr)           #
library(tidyr)           #
library(Hmisc)           #
##########################


##################################################################################################
#                                load bioanalyzer data                                           #
##################################################################################################
# relative number of fragments of input fragment distributions as measured by bioanalyzer        #
d.df<- fread(here::here("Marcelo_paper_figures_scripts_in_order",                                #
                        "supplementary_figures/S_fig1/bioanalyzer_frag_dist.csv"))               #
# pick only NP small fragment distribution                                                       #
up15<- d.df %>% dplyr::filter(sample == "UP15")                                                  #
# remove rows with NA values                                                                     #
up15.dist<- up15[1:695,]                                                                         #
##################################################################################################

###################################################################################################
# artificially extend the distribution to 50 since it is incomplete because the                   #
# bioanalyzer size control overlap the 50 - 70 bp mark                                            #
###################################################################################################
a <- seq(50,71.5, by = 0.3)                                                                       #
fi <- up15$rel_mol[1]                                                                             #
# make sure extended bases and evenly spread                                                      #
a0 <- (fi)/length(a)                                                                              #
a2 <- seq(0, fi, by = a0)                                                                         #
# add extended bases and relative fragments to the bioanalyzer dataset                            #
bases <- c(a, up15$bases[1:695])                                                                  #
rel_frag <- c(a2[1:72], up15$rel_mol[1:695])                                                      #
up15.dist <- data.table(bases = as.integer(bases), rel_frag = rel_frag)                            #
###################################################################################################


####################################################################################################
# fit a loess model line extracting 1 point for each base pair 
####################################################################################################
p1 <- ggplot(aes(x = bases, y = rel_frag), data = up15.dist) +
  geom_point() + 
  geom_smooth(method = "loess", span = 0.2, n = 569) + #span controls the amount of smoothing
  scale_x_continuous(limits = c(0,1200), breaks = c(0,100,300,500,700,1000,1200), expand = c(0, 0))+
  ggtitle("Fragment distribution of short input fragments") +
  labs(x = "base pairs", y = "Relative # of fragments") +
  theme_bw() +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text  = element_text(size=18),
        axis.title = element_text(size = 18, face = "bold"))
p1
# save a dataframe with 60 points if fitted loess y-hat data
data.d <- ggplot_build(p1)$data[[2]]
# make sure first point starts at 0
data.d$y[which(data.d$y < 0)] <- 0
####################################################################################################



####################################################################################################
#                     Calculate the frequency of a base pair fragment                              #
####################################################################################################
# calculate the frequency of the distribution based on the fitted relative number of fragments     #
density <- data.d$y/(sum(data.d$y))                                                                #
# check the sum is 1                                                                               #
sum(density) # should be 1                                                                         #
# add to a dataframe                                                                               #
up15.dist1 <- data.table(bases = data.d$x, density = density)                                      #
####################################################################################################

# Plot frequency distribution of base pairs
p1 <- ggplot(aes(x = bases, y = density), data = up15.dist1) +
  geom_point() +
  geom_line() + #span controls the amount of smoothing
  scale_x_continuous(limits = c(0,700), breaks = seq(0 , 700, 100), expand = c(0, 0))+
  scale_y_continuous(limits = c(0,0.005), expand = c(0, 0)) +
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

up15.dist_all<- up15.dist1[,-3]

# save the matrix of all fragment sizes 
folder.name = "./Marcelo_paper_figures_new_order/model_examination/adjust_dist/"
write.csv(up15.dist_all, paste(folder.name,"up15.dist_all_frag_sizes.csv", sep = ""))

####################################################################################################
#                             Calculate a matrix of frequency density classes                      #
####################################################################################################
# make a vector of size classes limits                                                             #
cuts <- seq(50,620, by = 10)                                                                       #
# group fragment sizes by size classes                                                             #
test <- cut2(up15.dist1$bases, cuts = cuts )                                                       #
up15.dist1$class_size <- test                                                                      #
# get the mean density per size class                                                              #
d <- up15.dist1 %>% group_by(class_size) %>% summarise(avg = mean(density))                        #
# check if the sum is around 1                                                                     #
sum(d$avg)                                                                                         #
# get the middle of each fragment size class                                                       #
size.class <- seq(55,615, by = 10)                                                                 #
# save in a dataframe of fragment size classes                                                     #
up15.dist_size.class <- data.frame(size.class = size.class, ave.density = d$avg)                   #
# save the matrix                                                                                  #
folder.name = "./Marcelo_paper_figures_new_order/model_examination/adjust_dist/"                   #
write.csv(up15.dist_size.class, paste(folder.name,"up15.dist_size.class.csv", sep = ""))           #
####################################################################################################

up15.dist_size.class <- fread(paste(folder.name,"up15.dist_size.class.csv", sep = ""))

# plot fragment size classes
p1 <- ggplot(aes(x = size.class, y = ave.density), data = up15.dist_size.class) +
  geom_bar(stat= "identity") +
  scale_x_continuous(limits = c(0,700), breaks = seq(0 , 700, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,0.005), expand = c(0, 0)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
      panel.grid.minor = element_line(colour="white", size=0.5),
      axis.text  = element_text(size=10, face = "bold"),
      axis.title = element_text(size = 16, face = "bold")) 
p1


