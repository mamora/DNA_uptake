###########################################################
######              modeling functions             ########         
###########################################################

# This script list 1 function  

# This calculates the contribution to uptake from one fragment size
# in turn with a given proportion of the fragment distribution taking into account probability of encounter a fragment,
# probability of binding a fragment and probability of taking up a fragment

# later this function will be used on a loop to calculate predicted uptake from several fragment sizes for the entire genome

###################################################################################
#                         Parameters used in all functions                        #
#  st is the start genome position                                                #
#  end is the end genome position                                                 #
#  a = genome position                                                            #
#  b = list of USS positions                                                      #
#  w = fragment size in bases                                                     #
#  fr = proportion of a given fragment size in the entire distribution.           #
#       (example: "1" means 100% of the fragments are from that given size "w").  #
#  sd = sliding distance                                                          #
###################################################################################



##########################################
#  calculations needed to run prob.model #
##########################################


#make vectors that indicate probability of binding "aa0" vs base pair size "aa1"
aa0<- seq(0 , 1, length.out = 11975)
aa1<- seq(1 , 11975, 1)
str(t)




###############################################################################################
# calculate a vector of adjusted probabilities "adj" based on distance between USS            #
                                                                                              #    
USS_separation<- seq(0, 11975, by = 1)  # Adjust this to the max range of fragments sizes     #
                                                                                              #    
adjustment_2<-  ((USS_separation/sd)*1)/max(USS_separation/sd) # for 2 USS/frag               #
adjustment_3<-  ((USS_separation/sd)*2)/max(USS_separation/sd) # for 3 USS/frag               #
adjustment_4<-  ((USS_separation/sd)*3)/max(USS_separation/sd) # for > 3 USS/frag             #
###############################################################################################




                                              ##################
                                              # model function #
###########################################################################################################################################################
# model that calculates predicted uptake for a given fragment size with a given proportion on the fragment distribution.                                  #
#     in  this model probability of encounter a fragment, binding a fragment and taking up a fragment                                                     #
#      are calculated. Predicted uptake is calculated by multipliying the 3 probabilities                                                                 #      
###########################################################################################################################################################
# all the following code belong to the model function

prob.model_bugs_fixed<- function(a, b = circle.uss.list$keypos, w = 200, fr = 1, sd = 100){

# a is the focal position  
# b is a circularized list of USS10
# w is fragment size  
# fr is contribution to uptake of a given fragment size  
# sd is sliding distance    
  
################################################################################  
# define all positions included within the fragment size w from the position a #
  a1<- c((a):(a + (w - 1)))                                                    #
  a2<- c((a - 1):(a - (w - 1)))                                                #
  # which USSs are within a1 and a2  
  c.t<- b[b %in% a1]                                                           #
  c.b<- b[b %in% a2]                                                           #
  ctb<- c(c.t, c.b)                                                            #
  b1<- (ctb - a)  # b1 has the position with USS fragments                     #                                         #
  ##############################################################################
  
  
  
  #######################################################################################
  # calculate the number of fragments that will contain a full USS in the entire range  #   
   b2<- w - (abs(b1) + 15)                                                              # 
                                                                                        #
  # make a matrix of number of fragments with a USS when USS is close to the focal pos  #
  # when the central USS position of that USS is 16 bases or closer                     #
  s<- c(-16:16)                                                                         #
  v1<- w- c(31, rep(30, 31),31)                                                         #
  v3<- cbind(s, v1)                                                                     #
  v4<- b1[which(b1 %in% v3[,1])]                                                        #
  v5<- as.numeric(v3[which(v3[,1] %in% v4),2])                                          #
  b2[which(b1 %in% v3[,1])] <- v5                                                       #
  #######################################################################################
  
  
###########################################################################  
  # if there is no USS then assign background level predicted uptake      #
  if(length(ctb) == 0) {                                                  #
    total<- (0.003 * w) * (0.01) * fr # Bo = (0.003 * w), Uo = (0.01)     #
    return(total)                                                         #
###########################################################################  
    
    } else {

  # If there are USSs then calculate probability of encounter, binding and of being take up.

      
##################################################################################################################################      
########                        if there is is only 1 USS in the entire range                                               ######
                                                                                                                                 #  
      if(length(ctb) == 1) {    # (Buss * Uss + B0 * U0) * fr                                                                    #                            
                                                                                                                                 #        
        if(sd >= w){  # if sliding distance is > than fragment size then the binding protein should find the USS always          #
        total<- fr *((b2/w * circle.uss.list$p.taken.up[circle.uss.list$keypos %in% ctb] * 1) +   ((w - b2)* (3e-5 * w))/w)      # 
        # B0 * U0 = (3e-5 * w), number of fragments without a uss = (w - b2), b2 = fragments with a uss                          #
        # Buss =  b2/w, fr = contribution to uptake                                                                              #                            
        } else {   # if sliding distance is > than fragment size then binding prob depend on difference between sd/w             #
        total<- fr * (((b2 * circle.uss.list$p.taken.up[circle.uss.list$keypos %in% ctb] * sd/w) + ((w - b2)* (3e-5* w))/w))/w   #  
        }                                                                                                                        #  
        return(total)                                                                                                            #
      }                                                                                                                          #
##################################################################################################################################
    
    
    
##################################################################################################################################
#                                                     If there is more than 1 USS in the fragment range                          #               
                                                                                                                                 #        
# if there is more than 1 USS, calculation are more complicated and probabilities have to be estimatedfor each unique fragment   #
# alignment                                                                                                                      #  
                                                                                                                                 #        
    # make a matrix with all unique fragment starting and end positions                                                          #  
    s<- (a - (w - 1))                                                                                                            #
    m<- matrix(nrow = w, ncol =1)                                                                                                #
    m[1,1]<- s                                                                                                                   #
                    
    # add first position of each unique fragment alignment                                                                       #        
    for(i in 1:(w - 1)){                                                                                                         #
      m[(i + 1),1]<- s + i                                                                                                       #
    }                                                                                                                            #
                                                                                                                                 #  
    num.uss<- rep(0, w) # add an empty column                                                                                    #
    
    # which USS are present in each unique fragment alignment                                                                    #  
    for(i in 1:w){                                                                                                               #
      temp<-  ctb[ctb >= (m[i,1] + 14) & ctb <= ((m[i,1] + (w - 1)) - 14)]                                                       #
      num.uss[i]<- length(temp)                                                                                                  #
    }                                                                                                                            #
    one.or.more<- which(num.uss >= 1)                                                                                            #
                                                                                                                                 #  
                                                                                                                                 #  
    uss<- rep(0, w) # add an empty column                                                                                        #  
    p.bind<- rep((0.003 * w), w) # add an column with B0 , will be relaced if there is a USS                                     #  
    separ<- rep(0, w) # add an empty column                                                                                      #
    p.taken.up<- rep(0.01, w) # add an column with U0 , will be relaced if there is a USS                                        #      
                                                                                                                                 #
    for(i in one.or.more){    # calculate Buss and Uuss for each fragment alignment                                              #                                                                                                      #
      temp<-  ctb[ctb >= (m[i,1] + 14) & ctb <= ((m[i,1] + (w - 1)) - 14)]                                                       #
                                                                                                                                 #    
      if(sd >= w){                                                                                                               #
        p.bind[i]<- 1                                                                                                            #
      } else{                                                                                                                    #
        p.bind[i]<- sd/w                                                                                                         #
      }                                                                                                                          #
                                                                                                                                 #
      if(length(temp) == 1){  # for fragment alignments with 1 USS                                                               #        
        uss[i]<- circle.uss.list$USS.score[which(circle.uss.list$keypos %in% temp)]                                              #
        p.taken.up[i]<- d1[which(abs(j- uss[i])==min(abs(j- uss[i])))]                                                           #
        
      } else {                                                                                                                   #
        if(length(temp) >= 2){  # if there is more than 1 uss in a fragment alignment                                            #
          # calculate separation between more distant USS in a fragment alignment cluster                                        #
          separ<- abs(min(temp) - max(temp))                                                                                     #
          #if 2 or more USS are found in the unique fragment calculate the mean USS score                                        #
          uss[i]<- mean(c(circle.uss.list$USS.score[circle.uss.list$keypos %in% temp]))                                          #
          p.taken.up[i]<- d1[which(abs(j- uss[i])==min(abs(j- uss[i])))]                                                         #
          # add adjustment, based on uss separation, to probability of binding to the fragment                                   #  
          # adjustment value added is based n the number of USS                                                                  #
          if (length(temp) == 2) {                                                                                               #
            p.bind[i]<- p.bind[i] + adjustment_2[USS_separation %in% separ]                                                      #  
                                                                                                                                 #
          } else {                                                                                                               #
            if (length(temp) == 3) {                                                                                             #
              p.bind[i]<- p.bind[i] + adjustment_3[USS_separation %in% separ]                                                    #
                                                                                                                                 #                                                        
            } else {                                                                                                             #   
              if (length(temp) >= 4) {                                                                                           #
                p.bind[i]<- p.bind[i] + adjustment_4[USS_separation %in% separ]                                                  #
              }                                                                                                                  #
            }                                                                                                                    #                                                                                                                    #  
          }                                                                                                                      #
        }                                                                                                                        #
      }                                                                                                                          #
    }                                                                                                                            #      
     # sum Uuss and Buss of all fragment alignments and multiply that for fragment size contribution                             #                                                                                                
    total<-  (sum(p.bind * p.taken.up)/w) * fr                                                                                   #  
##################################################################################################################################                
    
    return(total) 
    
  }
  
} 

# compilate in C++ the function to increase speed
prob.model_f <- compiler::cmpfun(prob.model_bugs_fixed)