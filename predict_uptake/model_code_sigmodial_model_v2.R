###########################################################
######              modeling functions             ########         
###########################################################

#############################################################################
#                         Parameters used in all functions                  #
#  a = genome position                                                      #
#  b = list of USS positions                                                #
#  w = fragment size in bases                                               #
#  b0 = binding relative frequency when a fragment lacks a uss              #
#  u0 = uptake initiation when a fragment lacks a uss                       #
#############################################################################


#####################################################################################
#   calculate expected uptake given sigmodial curve of observed uptake vs USS score #
#####################################################################################
# Imax =  maximum intensity 3.53                                                    #
Imax = 3.53                                                                         #
# a1 = slope at tmid 3.62                                                           #
a1 = 3.62                                                                           #
# tmid =  time at half intensity  10.6                                              #
tmid = 10.6                                                                         #
# t = time (x parameter)                                                            #
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5774301/                             #
# time is equal  to the USS_score                                                   #
# d1 is equal to the expected uptake                                                #
time <- seq(9, 13, by = 0.01)                                                       #
d1<- c()                                                                            #
for (i in 1:length(time)){                                                          #
  d0<- Imax/(1 + exp(-a1 *(time[i] - tmid)))                                        #
  d1<- c(d1, d0)                                                                    #
}                                                                                   #
#####################################################################################

##################
# model function #
################################################################################
# This function generates a matrix of binding and uptake initiation            #
# probabilities for all posible fragment aligments containing a focal position #
################################################################################
# the following code belong to the model function up to the "end" comment
prob.model_sig_v2 <- function(a, b = b, w = w, b0 = b0, u0 = u0 ){
 
  # Calculate the position of USS present in all possible fragment alignment
  # containing focal position "a"
  
  #### a1 = extreme left end of region under consideration
  #### a2 = extreme right end of region under consideration
  #### ctb = list of USS10 positions in the region under consideration
  
  a1 <- c((a):(a + (w - 1)))                                                    
  a2 <- c((a - 1):(a - (w - 1)))                                                
  c.t <- b[b %in% a1]                                                           
  c.b <- b[b %in% a2]                                                           
  ctb <- c(c.t, c.b)
  
  # Make a matrix with the starting position of each fragment alignment
  # containing focal position "a"
  s <- (a - (w - 1))                                                            
  m <- matrix(nrow = w, ncol = 1)                                               
  m[1,1] <- s
  for(i in 1:(w - 1)){                                                          
    m[(i + 1),1] <- s + i                                                       
  }
  
  # Calculate the number of USSs present in each unique fragment alignment      

      num.uss <- rep(0, w)   # Create the empty vector 'num.uss'                
      # loop calculating number of USS10s in each fragment of matrix m
      # temp holds the USS10 positions for each fragment in turn
      # The '14' calculations exclude partial USSs at the fragment ends         
      for(i in 1:w){           
      temp<-  ctb[ctb >= (m[i,1] + 14) & ctb <= ((m[i,1] + (w - 1)) - 14)]      
      num.uss[i] <- length(temp)                                                
    }                                                                           
    
      # list of rows in num.uss and m that have 2 or more USS10s
    one.or.more <- which(num.uss >= 1)                                          
    
    # Create 4 vectors to contain (1) USS10 score (2) binding probability, 
    # (3) USS10 separation and (4) uptake initiation probability                 
    uss <- rep(0, w)           # create an empty vector for USS10 scores        
    p.bind <- rep((b0), w)     # add a vector with baseline binding probability
    # it will be replaced if there is a USS                                    
    separ <- rep(0, w)         # create an empty vector for USS10 separation                       
    p.taken.up <- rep(u0, w)   # add a vector with baseline uptake probability
    # it will be replaced if there is a USS
    
    # For each fragment alignment, this loop calculate binding 
    # and uptake probability
    for(i in one.or.more){                                                      
      temp <-  ctb[ctb >= (m[i,1] + 14) & ctb <= ((m[i,1] + (w - 1)) - 14)]      
      if(length(temp) == 1){  # for fragment alignments with 1 USS              
        p.bind[i] <- 1 - w/28000  
        uss<- circle.uss.list$USS.score[which(circle.uss.list$keypos %in% temp)]
        p.taken.up[i]<- d1[which(time %in% uss)]                                
        } else {                                                                  
        if(length(temp) >= 2){ # for fragment alignments with 2 or more USS    
          separ <- abs(min(temp) - max(temp))  # distance of the 2 farther uss    
          uss <- (circle.uss.list$USS.score[circle.uss.list$keypos %in% temp])   
          # calculate frequency of uptake 
          p.taken.up[i]<- mean(d1[which(time %in% uss)])              
          p.bind[i] <- 1 - ((w/28000)*(1- (separ/w)))                          
        
        }                                                                  
      }                                                                       
    }                                                                         
    # join all vectors of binding and uptake probabilities into a matrix "m1" 
    m1 <- cbind(m, p.bind,p.taken.up)
    exp <- (sum(m1[,2] * m1[,3]))
    L <- exp - (sum(m1[1,2] * m1[1,3]))
    if(L < 0){
    L <- 0
    }
    newList <- list("exp" = exp, "L" = L)
return(newList)    
}                                                                             
# end of the function


left_pe <- function(a, b = b, w = w, b0 = b0, u0 = u0, newList){  
  
  # Add a fragment alignment to the focal position next to the one calculated 
  # in matrix "m1" and see if it has USS
  a1 <- a 
  a2 <- c((a1 - w + 1)):(a1)                                                    
  ctb <- b[b %in% a2]                                                           
  
  # temp holds the USS10 positions for the fragment alignment
  # The '14' calculations exclude partial USSs at the fragment ends          
  temp <- ctb[ctb >= (a2[1] + 14) & ctb <= ((a2[1] + (w - 1)) - 14)]
  
  # calculates binding and uptake probability if the fragment has no USS   
  if (length(temp) == 0) {
    p.b <- b0
    p.t <- u0
    
    # calculates binding and uptake probability if the fragment has 1 USS
  } else if (length(temp) == 1) {
    p.b <- 1 - (w/28000)                                                     
    uss <- circle.uss.list$USS.score[which(circle.uss.list$keypos %in% temp)] 
    p.t <- d1[which(time %in% uss)] 
    
    # calculates binding and uptake probability if the fragment has 2 or more USSs
  } else {
    separ <- abs(min(temp) - max(temp))                                    
    uss <- (circle.uss.list$USS.score[circle.uss.list$keypos %in% temp])     
    p.t <- mean(d1[which(time %in% uss)])                                  
    p.b <- 1 - ((w/28000)*(1- (separ/w)))  
  }
  
  Ln <- p.b * p.t
  return(Ln)
}



#################################               
#       loop function           #   
################################################################################
# This function takes the matrix from prob.model_simple_v1.1 and calculates    #
# binding and uptake prob for next sliding window focal position               #
################################################################################
# the following code belong to the model function up to the "end" comment
right_pe <- function(a, b = b, w = w, b0 = b0, u0 = u0, exp = newList$exp, L = newList$L){  

# Add a fragment alignment to the focal position next to the one calculated 
# in matrix "m1" and see if it has USS
  a1 <-  a
  a2 <- c((a1):(a1 + (w - 1)))                                                    
  ctb <- b[b %in% a2]                                                           

# temp holds the USS10 positions for the fragment alignment
# The '14' calculations exclude partial USSs at the fragment ends          
temp <- ctb[ctb >= (a1 + 14) & ctb <= ((a1 + (w - 1)) - 14)]

# calculates binding and uptake probability if the fragment has no USS   
  if (length(temp) == 0) {
      p.b <- b0
      p.t <- u0

# calculates binding and uptake probability if the fragment has 1 USS
          } else if (length(temp) == 1) {
       p.b <- 1 - (w/28000)                                                     
      uss <- circle.uss.list$USS.score[which(circle.uss.list$keypos %in% temp)] 
      p.t <- d1[which(time %in% uss)] 

# calculates binding and uptake probability if the fragment has 2 or more USSs
             } else {
        separ <- abs(min(temp) - max(temp))                                    
        uss <- (circle.uss.list$USS.score[circle.uss.list$keypos %in% temp])     
        p.t <- mean(d1[which(time %in% uss)])                                   
        p.b <- 1 - ((w/28000)*(1- (separ/w)))  
        }

Rn <- p.b * p.t
new.exp<- L + Rn
Ln <- left_pe(a = a, b = b, w = w, b0 = b0, u0 = u0)
new.L <- new.exp - Ln
if(new.L < 0){
  new.L <- 0
}
newList <- list("exp" = new.exp, "L" = new.L)
return(newList)
} 
# end of the function

################################
# Matrix generation function   #
################################################################################
# function to calculate a list of tables where each table has binding and      #
# uptake probabilities of a set of genome positions in the vector "genome"     #
################################################################################
func.model <- function(ds){

  sim <- vector() # create empty vector

  # generate a matrix of uptake and binding probabilities 
  # for first focal position
  newList<- prob.model_sig_v2(a = genome[1], b = b, w = up15.dist1$bases[ds],
                              b0 = b0, u0 = u0)
  sim<- c(sim, newList$exp)
 
  # use the matrix containing binding and uptake probabilities from the 
  # previous focal position and calculate binding and uptake probability 
  # of the next focal position.
  for (i in 2:length(genome)){
    newList <- right_pe(a = genome[i], b = b, w = up15.dist1$bases[ds], b0 = b0, u0 = u0,
                        exp = newList$exp, L = newList$L )
    sim<- c(sim, newList$exp)
  }
  return(sim)
}





