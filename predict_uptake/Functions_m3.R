
###########################################################
##     Modeling functions used by DNA uptake model       ##         
###########################################################


##############################################################################
# HOW THE UPTAKE PREDICTIONS WORK:
# The model does a full calculation of uptake contributions all the overlapping 
# fragments only for the first focal position (in function prob.model_sigv2).
# For each subsequent focal position, the model does not calculate the contributions 
# of all the overlapping fragments. Instead it take the previous position's sum and 
# modifies it by adding the contribution of the new rightmost fragment
# and subtracting the contribution of the former leftmost fragment (in function right_pe).
# It cannot efficiently save the contribution of the formerly leftmost fragment so 
# it recalculates it (in function left_pe). Function model first calls prob.model_sigv2 
# to calculate the first position, and then iterates through the other positions using function right_pe.
##############################################################################


#############################################################################
#                         Parameters used in all functions                  #
#  focal = Focal genome position                                            #
#  USS_list = list of USS positions                                         #
#  L = fragment length in base pairs                                        #
#  b0 = binding relative frequency when a fragment lacks a uss              #
#  u0 = uptake initiation when a fragment lacks a uss                       #
#############################################################################

##################################
##  Function uptake_first_pos   ##
################################################################################
# Given a specified size-class of DNA fragment, this function generates the    #
# sum of the uptake contributions for DNA fragments of this size that overlap  #
# the first focal position in the DNA segment or genome being analyzed.        #
# Subsequent functions iteratively generate the contribution sums for          #
# all the other positions in the segment or genome, for the same size class.   #
################################################################################

uptake_first_pos <- function(focal, USS_list = USS_list, L = L, b0 = Baseline_binding, u0 = Baseline_uptake ){
  
  # Specify the ends of the sub-segment containing all the this-size fragments 
  # that overlap the focal position "focal"
  
  #### a1 = extreme left end of region under consideration
  #### a2 = extreme right end of region under consideration
  #### list = list of USS positions in the region under consideration  
  
  a1 <- c((focal):(focal + (L - 1)))                                                    
  a2 <- c((focal - 1):(focal - (L - 1)))                                                
  list_right <- USS_list$pos[USS_list$pos %in% a1]                                                           
  list_left <- USS_list$pos[USS_list$pos %in% a2]                                                           
  list <- c(list_right, list_left)
  
  # Make a matrix 'm' with the starting position of each overlapping fragment 
  
  s <- (focal - (L - 1))                                                            
  m <- matrix(nrow = L, ncol = 1)                                               
  m[1,1] <- s
  for (i in 1:(L - 1)) {                                                          
    m[(i + 1),1] <- s + i                                                       
  }
  
  # Calculate the number of USSs present in each overlapping fragment.
  
  # Create the empty vector 'num.uss'
  num.uss <- rep(0, L)     
  
  # loop calculating the number of USSs in each fragment of matrix m.
  # temp holds the USS positions for each fragment in turn.
  # The '14' calculations exclude partial USSs at the fragment ends.         
  for (i in 1:L) {           
    temp <-  list[list >= (m[i,1] + 14) & list <= ((m[i,1] + (L - 1)) - 14)]      
    num.uss[i] <- length(temp)                                                
  }                                                                           
  
  # list of rows in num.uss and m that have 1 or more USSs
  one.or.more <- which(num.uss >= 1)                                          
  
  # Create 4 vectors to contain 
  # (1) USS score ('uss'), 
  # (2) binding probability ('p_bind'), 
  # (3) USS separation (if more than 1) ('separ')
  # (4) uptake initiation probability ('p_taken_up')     
  
  uss <- rep(0, L)           # Create an empty vector for USS scores.        
  p_bind <- rep((b0), L)     # Add a vector with baseline binding probability.
                               # It will be replaced if there is a USS                                    
  p_taken_up <- rep(u0, L)   # Add a vector with baseline uptake probability.
                                # It will be replaced if there is a USS.
  
  ###############################################################################################
  ##  This is the heart of the model.  Here we specify how to calculate the specific-          ##  
  ##  binding and uptake probabilities for each fragment as a function of their USSs.          ##
  ###############################################################################################
  
  # Overlapping fragments that have no USS have been assigned the baseline 
  # binding (p_bind) and uptake (p_taken_up) probabilities.
  # For each overlapping fragment that has one USS, this loop first assigns a p_bind value
  # (high for short fragments, lowered for long fragments in proportion to the maximimum fragment length).
  # It then extracts the relevant USS score from the terminally-redundant USS list and used this to get
  # p_taken_up from the Uptake_table.  
  # For fragments with 2 or more USS, it first calculates the separation between the 
  # most-distant USSs.  It then calculates p_bind, lowered for long fragments in proportion to both the 
  # length of the fragment and the closeness of the USSs.  It then extracts the relevant USS score from the 
  # terminally-redundant USS list, gets the corresponding uptake probabilities values
  # from the Uptake_table. and calculates p_taken_up as the max of thes  uptake probabilities.

  for (i in one.or.more) {                                                      
    temp <-  list[list >= (m[i,1] + 14) & list <= ((m[i,1] + (L - 1)) - 14)]      
 # for fragments with 1 USS 
       if (length(temp) == 1) {                
      p_bind[i] <- p_b_uss$p_b[1] 
      uss<- USS_list$USS.score[which(USS_list$pos %in% temp)]
      p_taken_up[i] <- Uptake_table$predicted_uptake[which(Uptake_table$score %in% uss)]                                
    } else {  
 # for fragments with 2 or more USS 
      if (length(temp) >= 2) {    
        p_bind[i] <- p_b_uss$p_b[length(temp)]         
        uss <- USS_list$USS.score[USS_list$pos %in% temp]   
        # calculate frequency of uptake 
        p_taken_up[i] <- mean(Uptake_table$predicted_uptake[which(Uptake_table$score %in% uss)])              
      }                                                                  
    }                                                                       
  }        
  
  #  At this poinbt we have calculated the binding and uptake predictions for 
  #  every fragment overlapping the focal position.   
  
  
  # join all vectors of binding and uptake probabilities into a matrix "m1" 
  m1 <- cbind(m, p_bind,p_taken_up)
  
  # For each fragment, calculate the total uptake as the product of p_bind and p_taken_up, 
  # and sum these to get 'exp', the total uptake for this position from this size class of fragments.
  
  exp <- (sum(m1[,2] * m1[,3]))
  
  # Prepare to calculate the new 'exp' for the next focal position.
  # Subtract from the current 'exp' the total uptake due to the leftmost fragment, 
  # since this fragment will not contribute to uptake at the next focal position.
  # Save this value as 'overlap_frag'.
  overlap_frag <- exp - (sum(m1[1,2] * m1[1,3]))
    if (overlap_frag < 0) {
      overlap_frag <- 0
    }
  
# Create 'newlist' which contains both 'exp and 'overlap_frag'
# Return 'newlist' to function 'model', which called it.
  newList <- list("exp" = exp, "overlap_frag" = overlap_frag)
  return(newList)    
} 

#  End of Function prob.model_sigv2
##############################################################################



#################################               
#    Function  left_pe       
################################################################################
# For each new focal position, this function recalculates the contribution of the leftmost fragment
# so this value can be subtracted from the next position's contribution.
##############################################################################

left_pe <- function(focal, USS_list = USS_list, L = L, b0 = Baseline_binding, u0 = Baseline_uptake, newList){  
  
  # Add a fragment alignment to the focal position next to the one calculated 
  # in matrix "m1" and see if it has USS
  a1 <- focal 
  a2 <- c((a1 - L + 1)):(a1)                                                    
  list <- USS_list$pos[USS_list$pos %in% a2]                                                           
  
  # temp holds the USS10 positions for the fragment alignment
  # The '14' calculations exclude partial USSs at the fragment ends          
  temp <- list[list >= (a2[1] + 14) & list <= ((a2[1] + (L - 1)) - 14)]
  
  # calculates binding and uptake probability if the fragment has no USS   
  if (length(temp) == 0) {
    p_bind <- b0
    p_taken_up <- u0
    
    # calculates binding and uptake probability if the fragment has 1 USS
  } else if (length(temp) == 1) {
    p_bind <- p_b_uss$p_b[1]
    uss <- USS_list$USS.score[which(USS_list$pos %in% temp)] 
    p_taken_up <- Uptake_table$predicted_uptake[which(Uptake_table$score %in% uss)] 
    
    # calculates binding and uptake probability if the fragment has 2 or more USSs
  } else {
    uss <- USS_list$USS.score[USS_list$pos %in% temp]     
    p_taken_up <- mean(Uptake_table$predicted_uptake[which(Uptake_table$score %in% uss)])
    p_bind <- p_b_uss$p_b[length(temp)] 
  }
  
  # For the leftmost fragment, calculate the total uptake as the product of p_bind and p_taken_up, 
  Ln <- p_bind * p_taken_up
  return(Ln)
}

#  End of Function left_pe
##############################################################################


#################################               
#    Function  right_pe   #   
################################################################################
# This function calculates the uptake contribution of the new rightmost fragment for each new focal position, and adds it
# to the 'overlap_frag' value from the previous position (sum of all contributions except the formerly leftmost 
# fragment) to get the 'new_exp', the total contribution for this position (for this size class).
# It then calls left_pe to calculate the uptake contribution of the current leftmost fragment, so this value 
# can be subtracted from new_exp to generate the 'new_L' for the next position. for next sliding window focal position             
################################################################################

right_pe <- function(focal, USS_list = USS_list, L = L, b0 = Baseline_binding, u0 = Baseline_uptake, exp = newList$exp, overlap_frag = newList$overlap_frag){  
  
  # Add a rightmost-overlapping fragment  to the focal position and see if it has USS
  a1 <-  focal
  a2 <- c((a1):(a1 + (L - 1)))                                                    
  list <- USS_list$pos[USS_list$pos %in% a2]                                                           
  
  # temp holds the USS10 positions for the fragment alignment
  # The '14' calculations exclude partial USSs at the fragment ends          
  temp <- list[list >= (a1 + 14) & list <= ((a1 + (L - 1)) - 14)]
  
  # calculates binding and uptake probability if the fragment has no USS   
  if (length(temp) == 0) {
    p_bind <- b0
    p_taken_up <- u0
    
    # calculates binding and uptake probability if the fragment has 1 USS
  } else if (length(temp) == 1) {
    p_bind <- p_b_uss$p_b[1]                                                     
    uss <- USS_list$USS.score[which(USS_list$pos %in% temp)] 
    p_taken_up <- Uptake_table$predicted_uptake[which(Uptake_table$score %in% uss)] 
    
    # calculates binding and uptake probability if the fragment has 2 or more USSs
  } else {
    uss <- USS_list$USS.score[USS_list$pos %in% temp]     
    p_taken_up <- mean(Uptake_table$predicted_uptake[which(Uptake_table$score %in% uss)])
    p_bind <- p_b_uss$p_b[length(temp)] 
  }
  
  # Adds the new fragment's uptake contribution (prodict of binding and uptake)
  # to the previously summed contribution of all the other fragments ('overlap_frag' from the previous cycle)to get new_exp.
  Rn <- p_bind * p_taken_up
  new.exp <- overlap_frag + Rn
  
  # Now it calls function left_pe to recalculate the contribution of the leftmost fragment
  # so that can be subtracted from new_exp to get new_L.
  Ln <- left_pe(focal = focal, USS_list = USS_list, L = L, b0 = Baseline_binding, u0 = Baseline_uptake)
  new.overlap_frag <- new.exp - Ln
  if (new.overlap_frag < 0) {
    new.overlap_frag <- 0
  }
  
  # It returns new_list which provides exp and overlap_frag for the next position.
  newList <- list("exp" = new.exp, "overlap_frag" = new.overlap_frag)
  return(newList)
} 
# end of Function right_pe
################################################################################


################################
#Function all_genome_uptake  #
################################################################################
# function to calculate a vector 'sim  containing the contributions of fragments of the 
# specified size class to all positions in the genome or segment.
################################################################################
all_genome_uptake <- function(ds){
  
  sim <- vector() # create empty vector
  
  # Call function uptake_first_pos to calculate the total contribution to uptake of all   
  # the fragments of the specified size class that overlap the first focal position
  newList <- uptake_first_pos(focal = genome[1], USS_list = USS_list, L = Fragment_distribution_file$bases[ds],
                              b0 = Baseline_binding, u0 = Baseline_uptake)
  
  # Adds the contribution for the first position to a vector 'sim' that will ahve the spoistion-specific uptake values.
  sim <- c(sim, newList$exp)
  
  # Call function right_pe to calculate the total contribution to uptake of all   
  # the fragments of the specified size class that overlap each subsequent focal position.
  # right_pe will also call left_pe.)
  for (i in 2:length(genome)){
    newList <- right_pe(focal = genome[i], USS_list = USS_list, L = Fragment_distribution_file$bases[ds], b0 = Baseline_binding, u0 = Baseline_uptake,
                        exp = newList$exp, overlap_frag = newList$overlap_frag )
    
    # Adds each new position's contribution to the vector of position-specific uptakes.
    sim <- c(sim, newList$exp)
  }
  
  # Returns the vector 'sim' containing the contributions of fragments of the 
  # specified size class to all positions in the genome or segment.
  return(sim)
}
# end of Function all_genome_uptake
################################################################################