######################################################################
## create a model that predicts DNA uptake in NPlarge frag data #####
######################################################################

# the objective of this script is to run the model to predict uptake


####################################################################################
######             1. load samples and working directory                           #         
library(dplyr)                                                                     #
library(tidyr)                                                                     #
library(data.table)                                                                #                    #
                                                                                   #
# uptake ratios
Uptake.ratio.np<- fread("./datasets/final_datasets/DNA_uptake/Uptake.ratio.np.corrected.csv") #read uptake file 

# list of USS10
Up.USS.np.10.list<- fread(here::here("datasets/final_datasets","Uptake.uss.10.list.np.csv")) #save file

up15.dist1<- fread("./Marcelo_paper_figures_new_order/model/Sept_10_2018/up15.dist_adj_4.w.csv")


                                                                                    #
####################################################################################


# Settings

settings<- read.csv("./Marcelo_paper_figures_new_order/model/Sept_10_2018/Settings.csv")


###################################################
#                 Settings                        #
###################################################
# set the genomic region to be predicted          #
sd = settings$sd
b0 = settings$b0
u0 = settings$u0                                  #
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






#####################################################################################################################
#       incorporate the model function to a loop that calculates contribution to uptake                             #
#                               of all fragment sizes for all genomic positions                                     #  


b = circle.uss.list$keypos

source("./Marcelo_paper_figures_new_order/model/Sept_10_2018/model_code_v4.1.R")                                                                                       #


up15.dist1$bases<- as.integer(up15.dist1$bases)

func.model<- function(ds){
sim<- vector()
tab<- list()
ti<- prob.model_v4.1_p1(a = genome[1], b = b, w = up15.dist1$bases[ds] , sd = sd, b0 = b0, u0 = u0)
sum1<-  (sum(ti[,2] * ti[,3])) * up15.dist1$density[ds]
sim<- c(sim,sum1)
tab[[1]]<- ti
for(i in 2:length(genome)){
  t2<- prob.model_f_v4.1(m1 = tab[[i-1]], b = b, w = up15.dist1$bases[ds], sd = sd, b0 = b0, u0 = u0)  
  tab[[i]]<- t2
  }
return(tab)
}

sum.total<- function(tab, sim){
  sim<- (sum(tab[,2] * tab[,3])) * up15.dist1$density[ds]
  return(sim)
}


st = 1                                                                                                                              
end = 100000

genome<- c(st:end)

f.tem<- data.frame(matrix(nrow = length(genome), ncol = length(up15.dist1$bases)))

for(ds in 1:length(up15.dist1$bases)){ 
tab.1<- func.model(ds)   
sum<- sapply(X = tab.1, FUN = sum.total, sim = sim)
f.tem[,ds] <- sum 
}

write.csv(f.tem, file = "./Marcelo_paper_figures_new_order/model/Sept_10_2018/new_model_100kb_np_dist_w4.csv") 




                                                                                                                                                                                                                                      #
#####################################################################################################################