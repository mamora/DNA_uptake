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


###############################################################################################
# calculate a vector of adjusted probabilities "adj" based on distance between USS            #
#    
USS_separation<- seq(0, 600, by = 1)  # Adjust this to the max range of fragments sizes     #
#    
adjustment_2<-  ((USS_separation/sd)*1)/max(USS_separation/sd) # for 2 USS/frag               #
adjustment_3<-  ((USS_separation/sd)*2)/max(USS_separation/sd) # for 3 USS/frag               #
adjustment_4<-  ((USS_separation/sd)*3)/max(USS_separation/sd) # for > 3 USS/frag             #
###############################################################################################



##########################################################################################################
#                         Use sigmodial equation to calculate probability to be taken up                 #  
#

j<- seq(0, 12.7, by = 0.01) # sequence of USS scores                                                     #
Asym = 3.77 
xmid = 10.6
scal =0.15

d1<- c()
for (i in 1:length(j)){
  d0<- Asym/(1+exp((xmid-j[i])/scal))
  d1<- c(d1, d0) 
}                                                                                                        #


#adjust to  a range from 0 to 1                                                                          #
d1 <-  d1/Asym                                                                                        #

m<- data.frame(cbind(j, as.numeric(format(round(d1, 2), nsmall = 2))))


m.s<- m[(which(m$j < 11)),]

colnames(m.s)<- c("score", "uptake")

j2<- seq(11, 12.7, by = 0.01) # sequence of USS scores                                                     #

intercept = -1.1
slope = 0.4  

dl<- c()
for (i in 1:length(j2)){
  d0<- intercept + (slope *j2[i])
  dl<- c(dl, d0) 
}

dln<-seq(0.94, (0.94 + (1 - dl[1]/max(dl))), length.out = length(dl))

m2<- data.frame(cbind(j2, as.numeric(format(round(dln, 2), nsmall = 2))))

colnames(m2)<- c("score", "uptake")
m3<- rbind(m.s, m2)

circle.uss.list$USS.score<- as.numeric(format(round(circle.uss.list$USS.score, 2), nsmall = 2)) 

#
p.taken.up<- c()                                                                                         #
for(i in 1:length(circle.uss.list$USS.score)){                                                           #
  pu<- m3$uptake[which( m3$score %in% circle.uss.list$USS.score[i] )]                                    #
  p.taken.up<- c(p.taken.up, pu[1])                                                                      #
}                                                                                                        #

                                                                                                         #
circle.uss.list$p.taken.up<- p.taken.up                                                                  #
##########################################################################################################




##################
# model function #
###########################################################################################################################################################
# model that calculates predicted uptake for a given fragment size with a given proportion on the fragment distribution.                                  #
#     in  this model probability of encounter a fragment, binding a fragment and taking up a fragment                                                     #
#      are calculated. Predicted uptake is calculated by multipliying the 3 probabilities                                                                 #      
###########################################################################################################################################################
# all the following code belong to the model function


prob.model_v3.0_p1<- function(a, b = b, w = w, sd = sd, b0 = b0, u0 = u0 ){
  # a is the focal position  
  # b is a circularized list of USS10
  # w is fragment size  
  # fr is contribution to uptake of a given fragment size  
  # sd is sliding distance    
  a1<- c((a):(a + (w - 1)))                                                    #
  a2<- c((a - 1):(a - (w - 1)))                                                #
  # which USSs are within a1 and a2  
  c.t<- b[b %in% a1]                                                           #
  c.b<- b[b %in% a2]                                                           #
  ctb<- c(c.t, c.b)
  
  s<- (a - (w - 1))                                                                                                            #
  m<- matrix(nrow = w, ncol =1)                                                                                                #
  m[1,1]<- s                                                                                                                   #
  
  # add first position of each unique fragment alignment                                                                       #        
  for(i in 1:(w - 1)){                                                                                                         #
    m[(i + 1),1]<- s + i                                                                                                       #
  }
  
                                                                                                                                  
    num.uss<- rep(0, w) # add an empty column                                                                                    #
    
    # which USS are present in each unique fragment alignment                                                                    #  
    for(i in 1:w){                                                                                                               #
      temp<-  ctb[ctb >= (m[i,1] + 14) & ctb <= ((m[i,1] + (w - 1)) - 14)]                                                       #
      num.uss[i]<- length(temp)                                                                                                  #
    }                                                                                                                            #
    one.or.more<- which(num.uss >= 1)                                                                                            #
                                                                                                                                 #                                                                                                                                   #  
    uss<- rep(0, w) # add an empty column                                                                                        #  
    p.bind<- rep((b0), w) # add an column with B0 , will be replaced if there is a USS                                     #  
    separ<- rep(0, w) # add an empty column                                                                                      #
    p.taken.up<- rep(u0, w) # add an column with U0 , will be replaced if there is a USS                                        #      
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
        p.taken.up[i]<- d1[which(abs(j- uss[i])==min(abs(j- uss[i])))]                                                          #
                                                                                                                                 #
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
                p.bind[i]<- p.bind[i] + adjustment_4[USS_separation %in% separ]                                                #
              }                                                                                                                  #
            }                                                                                                                    #                                                                                                                    #  
          }                                                                                                                      #
        }                                                                                                                        #
      }                                                                                                                          #
    }                                                                                                                            #      

  m1<- cbind(m, p.bind,p.taken.up)
return(m1)    

}                                                                                           #  
    ##################################################################################################################################                
     
    
loop_the_genome<- function(m1, b = b, w = w, sd = sd, b0 = b0, u0 = u0){  

  a<- m1[w,1] + 1
  
  a1<- c((a):(a + (w - 1)))                                                    #
  a2<- c((a - 1):(a - (w - 1)))                                                #
  # which USSs are within a1 and a2  
  c.t<- b[b %in% a1]                                                           #
  c.b<- b[b %in% a2]                                                           #
  ctb<- c(c.t, c.b)
                                                                                                                    #   
temp<- ctb[ctb >= (a + 14) & ctb <= ((a + (w - 1)) - 14)]
    
  if(length(temp) == 0) {
      p.b<- b0
      p.t<- u0
    } else if (length(temp) == 1) {

            if(sd >= w){                                                                                                               #
            p.b<- 1                                                                                                            #
            } else{                                                                                                                    #
            p.b<- sd/w                                                                                                         #
          }
    
      scor<- circle.uss.list$USS.score[which(circle.uss.list$keypos %in% temp)]                                              #
      p.t<- d1[which(abs(j- scor)==min(abs(j- scor)))] 
    } else {
      if(sd >= w){                                                                                                               #
        p.b<- 1                                                                                                            #
      } else{                                                                                                                    #
        p.b<- sd/w                                                                                                         #
      }
        # calculate separation between more distant USS in a fragment alignment cluster                                        #
        separ<- abs(min(temp) - max(temp))                                                                                     #
        #if 2 or more USS are found in the unique fragment calculate the mean USS score                                        #
        p.tt<- mean(c(circle.uss.list$USS.score[circle.uss.list$keypos %in% temp]))                                          #
        p.t<- d1[which(abs(j- p.tt)==min(abs(j- p.tt)))]                                                         #
        # add adjustment, based on uss separation, to probability of binding to the fragment                                   #  
        # adjustment value added is based n the number of USS                                                                  #
        if (length(temp) == 2) {                                                                                               #
          p.b<- p.b + adjustment_2[USS_separation %in% separ]                                                     #  
          #
        } else {                                                                                                               #
          if (length(temp) == 3) {                                                                                             #
            p.b<- p.b  + adjustment_3[USS_separation %in% separ]                                                    #
            #                                                        
          } else {                                                                                                             #   
            if (length(temp) >= 4) {                                                                                           #
              p.b<- p.b  + adjustment_4[USS_separation %in% separ]                                                  #
            }
         }
      }
    }  

m2<- m1[-1,]

row<- c(a, p.b, p.t)

m3<- rbind(m2,row)

return(m3)
} 

    
# compilate in C++ the function to increase speed
prob.model_f_v3.0 <- compiler::cmpfun(loop_the_genome)


