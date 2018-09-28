###########################################################
######              modeling functions             ########         
###########################################################

#############################################################################
#                         Parameters used in all functions                  #
#  a = genome position                                                      #
#  b = list of USS positions                                                #
#  w = fragment size in bases                                               #
#  b0 = binding probability when a fragment lacks a uss                     #
#  u0 = uptake initiation probability when a fragment lacks a uss           #
#############################################################################


##################
# model function #
################################################################################
# This function generates a matrix of binding and uptake initiation            #
# probabilities for all posible fragment aligments containing a focal position #
################################################################################
# the following code belong to the model function up to the "end" comment
prob.model_simple_v1.1 <- function(a, b = b, w = w, b0 = b0, u0 = u0 ){
 
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
        p.taken.up[i] <-  (uss - 10)/(12.64 - 10)                                
        } else {                                                                  
        if(length(temp) >= 2){ # for fragment alignments with 2 or more USS    
          separ <- abs(min(temp) - max(temp))  # distance of the 2 farther uss    
          uss <- (circle.uss.list$USS.score[circle.uss.list$keypos %in% temp])   
          # probability of uptake = mean((uss.score - 10)/(max(uss.score) - 10))
          # max(uss.score) = 12.64
          p.taken.up[i] <- mean((uss - 10)/(12.64 - 10))              
          p.bind[i] <- 1 - ((w/28000)*(1- (separ/w)))                          
        
        }                                                                  
      }                                                                       
    }                                                                         
    # join all vectors of binding and uptake probabilities into a matrix "m1" 
    m1 <- cbind(m, p.bind,p.taken.up)
return(m1)    
}                                                                             
# end of the function

#################################               
#       loop function           #   
################################################################################
# This function takes the matrix from prob.model_simple_v1.1 and calculates    #
# binding and uptake prob for next sliding window focal position               #
################################################################################
# the following code belong to the model function up to the "end" comment
loop_the_genome_simple <- function(m1, b = b, w = w, b0 = b0, u0 = u0){  

# Add a fragment alignment to the focal position next to the one calculated 
# in matrix "m1" and see if it has USS
  a <- m1[w,1] + 1
  a1 <- c((a):(a + (w - 1)))                                                    
  ctb <- b[b %in% a1]                                                           

# temp holds the USS10 positions for the fragment alignment
# The '14' calculations exclude partial USSs at the fragment ends          
temp <- ctb[ctb >= (a + 14) & ctb <= ((a + (w - 1)) - 14)]

# calculates binding and uptake probability if the fragment has no USS   
  if (length(temp) == 0) {
      p.b <- b0
      p.t <- u0

# calculates binding and uptake probability if the fragment has 1 USS
          } else if (length(temp) == 1) {
       p.b <- 1 - (w/28000)                                                     
      uss <- circle.uss.list$USS.score[which(circle.uss.list$keypos %in% temp)] 
      p.t <- (uss - 10)/(12.64 - 10) 

# calculates binding and uptake probability if the fragment has 2 or more USSs
             } else {
        separ <- abs(min(temp) - max(temp))                                    
        uss <- (circle.uss.list$USS.score[circle.uss.list$keypos %in% temp])     
        p.t <- mean((uss - 10)/(12.64 - 10))                                   
        p.b <- 1 - ((w/28000)*(1- (separ/w)))  
        }

# delete the fragment aligment containing the previous focal position  
#that does not included the new focal position 
m2 <- m1[-1,]

# join all in a matrix 
row <- c(a, p.b, p.t)
m3 <- rbind(m2,row)

return(m3)
} 
# end of the function

########################################################################
#          compilate the function in C++ to increase speed             #
prob.model_f_simple_v1.1 <- compiler::cmpfun(loop_the_genome_simple)   #
########################################################################


################################
# Matrix generation function   #
################################################################################
# function to calculate a list of tables where each table has binding and      #
# uptake probabilities of a set of genome positions in the vector "genome"     #
################################################################################
func.model <- function(ds){

  sim <- vector() # create empty vector
  tab <- list() # create empty list

  # generate a matrix of uptake and binding probabilities 
  # for first focal position
  ti<- prob.model_simple_v1.1(a = genome[1], b = b, w = up15.dist1$bases[ds],
                              b0 = b0, u0 = u0)

    # Sum the product of binding per uptake probability of each fragment 
    # alignment. Multiply that per the frequency of the fragment size
  sum1 <-  (sum(ti[,2] * ti[,3])) * up15.dist1$density[ds]
  sim <- c(sim,sum1)
  tab[[1]] <- ti
  
  # use the matrix containing binding and uptake probabilities from the 
  # previous focal position and calculate binding and uptake probability 
  # of the next focal position.
  for (i in 2:length(genome)){
    t2 <- prob.model_f_simple_v1.1(m1 = tab[[i-1]], b = b,
                                  w = up15.dist1$bases[ds], b0 = b0, u0 = u0)  
    tab[[i]] <- t2
  }
  return(tab)
}

##########################################
# Predicted uptake calculation function  #
################################################################################
# For each table in the list "tab", Sum the product of binding per uptake      #
# probabilities and multiply that per fragment frequency to calculate          #
# predicted uptake                                                             #
################################################################################
sum.total <- function(tab, sim){
  sim <- (sum(tab[,2] * tab[,3])) * up15.dist1$density[ds]
  return(sim)
}



