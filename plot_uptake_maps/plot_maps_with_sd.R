###################################################################
#         plot uptake maps with standard deviation                #
###################################################################

######################################################
#                       WARNING                      #
#                                                    #
#   IN ORDER FOR THE SCRIPTS TO WORK CORRECTLY       #
#  YOU MUST FIRST OPEN THE PROJECT FILE NAMED:       #
#          " Uptake_summer2017.Rproj"                #
#       OTHERWISE EDIT THE "FOLDER.NAME" PATH TO     #
#   THE RIGHT LOCATION WHERE SCRIPTS ARE SAVED       #
######################################################

# LOAD PACKAGES
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)

# set scripts and output sigures folder name
folder.name <- ("./Marcelo_paper_figures_new_order/model/Models_plots/")

# set datasets folder name
folder.name2 <- ("./datasets/final_datasets/") 

# load USS scores
USS.scores<- fread(file = paste(folder.name2,
                                "USS_scores/uptake_model/np.USS.scores.corrected.csv", sep = ""))
#read uptake file
Uptake.ratio.np <- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv")      

# set USS arrows cutoff. If using sigmodial model then set arrow.uss.cut at 9.5
arrow.uss.cut <- 10

arrow_2.uss.cut <- 9.5

# remove extra columns added by fread function
Uptake.ratio.np$V1 <- NULL
Uptake.ratio.np$V1 <- NULL


########################################
#       set plots parameters           #
########################################


figure.title <- "Observed uptake ratios" # set figure tittle

x.axis.title <- "NP genomic positions" # set x axis tittle

y.axis.title <- "Uptake ratio" # set y axis tittle

y.axis.height<- 5  # set y axis top value. Not applicable for log 2 plot

y.axis.start<- 0  # set y axis top value. Not applicable for log 2 plot

start.pos <- 150000  # set starting position in the plot

end.pos <- 160000 # set end position in the plot
  
spacing_xaxis<- 2000  # sets spacing in x axis

arrows.y.height <- 0.2 # control the height of the USS arrows in the y axis

xy.title.size <- 22  # controls the title X and Y size

title.size <- 5  # controls the plot title size

legend.size<- 5 # controls the legend size

xy.num.text<- 24 # control xy numbers size

uptake<- Uptake.ratio.np$ratio_long # ratio of dataset to plot 

pos<- Uptake.ratio.np$pos # positions of the genome to plot

sd_short<- Uptake.ratio.np$sd_long

data<- data.table(pos, uptake, sd_short)

# load the functions to be used. MUST BE LOADED EVERY TIME A PARAMETER CHANGE
source(file = paste(folder.name,
                    "./plots_functions2.R",
                    sep = ""))

##################################################################################
#             plot standard uptake map. Predicted vs observed uptake             #
##################################################################################

# choose the theme to plot. 4 themes were included in the functions script:
#              theme_background_grey: includes a grey standard background
#              theme_background_white: includes a black and white background
#              theme_background_minimal: background without borders
#              theme_background_translucid: translucid background

# set add.arrows = TRUE to add USS arrows in the plot, otherwise set in FALSE

t<- standard.sd.plot (data = data, theme = theme_background_grey,  add.arrows = TRUE)

# ADD blue arrows for USS with scores 9.5 - 10. ignore this line if you do not want blue arrows

t2<- add_more_arrows(t)

output_name<- "peak_center24" # set plot figure name when its saved 

# save plot regular
file_name = paste(folder.name,output_name, ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, units = "px")
print(t2)
dev.off()

# save plot translucid
file_name = paste(folder.name,output_name, ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, bg = "transparent", units = "px")
print(t)
dev.off()

###############################################################################################
#             plot uptake map. Predicted vs observed uptake in scientific notation            #
###############################################################################################

# set add.arrows = TRUE to add USS arrows in the plot, otherwise set in FALSE

t<- scientif.plot (data = simple_model, theme = theme_background_grey,  add.arrows = TRUE)

t

output_name<- "test" # set plot figure name when its saved 

# save plot regular
file_name = paste(folder.name,output_name, ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, units = "px")
print(t)
dev.off()

# save plot translucid
file_name = paste(folder.name,output_name, ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, bg = "transparent", units = "px")
print(t)
dev.off()

#########################################################################################
#             plot log2 of y axis uptake maps. Predicted vs observed uptake             #
#########################################################################################

# log 2 of uptake ratios were plotted but legend in y axis shows the true uptake values

t<- log2.plot(data = simple_model, theme = theme_background_grey,  add.arrows = TRUE)

t

output_name<- "test" # set plot figure name when its saved 

# save plot regular
file_name = paste(folder.name,output_name, ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, units = "px")
print(t)
dev.off()

# save plot translucid
file_name = paste(folder.name,output_name, ".tiff", sep="")
tiff(file_name, width = 1200, height = 800, bg = "transparent", units = "px")
print(t)
dev.off()
