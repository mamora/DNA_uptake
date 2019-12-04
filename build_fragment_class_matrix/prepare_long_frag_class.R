#################################################################################################
#                   this script uses fragment distribution measured by bioanalyzer and          #
#                          then assign fragment sizes into 200bp classes                         #  
#################################################################################################

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
up13<- d.df %>% dplyr::filter(sample == "UP13")                                                  #
# remove rows with NA values                                                                     #
up13.dist<- up13[1:335,]                                                                         #
up13.dist <- data.table(bases = as.integer(up13.dist$bases), rel_frag = up13.dist$rel_mol)        #
###################################################################################################


####################################################################################################
# fit a loess model line extracting 1 point for each base pair 
####################################################################################################
seq <- c(846:14197)

p1 <- ggplot(aes(x = bases, y = rel_frag), data = up13.dist) +
  geom_point() + 
  geom_smooth(method = "loess", span = 0.2, n = 13352) + #span controls the amount of smoothing
  scale_x_continuous(limits = c(800,14500), expand = c(0, 0))+
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
up13.dist1 <- data.table(bases = data.d$x, density = density)                                      #
####################################################################################################

# Plot frequency distribution of base pairs
p1 <- ggplot(aes(x = bases, y = density), data = up13.dist1) +
  geom_point() +
  geom_line() + #span controls the amount of smoothing
  scale_x_continuous(limits = c(800,14500), breaks = seq(800 , 14500, 1000), expand = c(0, 0))+
  scale_y_continuous(limits = c(0,0.0004), expand = c(0, 0)) +
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

up13.dist_all<- up13.dist1

# save the matrix of all fragment sizes 
folder.name = "./Marcelo_paper_figures_new_order/model_examination/adjust_dist/"
write.csv(up13.dist_all, paste(folder.name,"up13.dist_all_frag_sizes.csv", sep = ""))

####################################################################################################
#                             Calculate a matrix of frequency density classes                      #
####################################################################################################
# make a vector of size classes limits                                                             #
cuts <- seq(800,14200, by = 200)                                                                   #
# group fragment sizes by size classes                                                             #
test <- cut2(up13.dist1$bases, cuts = cuts )                                                       #
up13.dist1$class_size <- test                                                                      #
# get the mean density per size class                                                              #
d <- up13.dist1 %>% group_by(class_size) %>% summarise(avg = mean(density))                        #
# check if the sum is around 1                                                                     #
sum(d$avg)                                                                                         #
# get the middle of each fragment size class                                                       #
size.class <- seq(900,14100, by = 200)                                                             #
# save in a dataframe of fragment size classes                                                     #
up13.dist_size.class <- data.frame(size.class = size.class, ave.density = d$avg)                   #
# save the matrix                                                                                  #
folder.name = "./Marcelo_paper_figures_new_order/model_examination/adjust_dist/"                   #
write.csv(up13.dist_size.class, paste(folder.name,"up13.dist_size.class.csv", sep = ""))           #
####################################################################################################

# plot fragment size classes
p1 <- ggplot(aes(x = size.class, y = ave.density), data = up13.dist_size.class) +
  geom_bar(stat= "identity") +
  scale_x_continuous(limits = c(800,14200), breaks = seq(800 , 14200, 2000), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,0.0005), expand = c(0, 0))
  p1


############################
# normalize to a sum of 1  #
############################
  
sum_all <-  sum(up13.dist_size.class$ave.density)
  
up13.dist_size.class$norm_freq <- up13.dist_size.class$ave.density/sum_all

sum(up13.dist_size.class$norm_freq)

NP_large_fragment_class <- data.table(bases = up13.dist_size.class$size.class,
                                      frequency = up13.dist_size.class$norm_freq)

write.csv(NP_large_fragment_class, paste(folder.name,"NP_large_fragment_class.csv", sep = ""))           #
