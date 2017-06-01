######################################################################
############## Functions for sites2pssm and scoring genomes ##########
######################################################################

###############
# Valid chars #
###############
chars <- c("a", "c", "g", "t", "n")
    # last character is reserved for ambiguities
    # this could be non-ACGTN, so amino acids, for example

chars.gap <- c("a", "c", "g", "t", "n","-")
# last character is reserved for ambiguities
# this could be non-ACGTN, so amino acids, for example


#############
# Utilities #
#############

# ggplot
library(ggplot2)

# Read in fastas, among other things
library(seqinr)

#sliding window
library(zoo)

# Split character string vectors into strings of characters
stringsplit <- function(string) unlist(strsplit(string, split=""))

# Make a table of counts from a string for membership in chars
tablepos <- function(string) table(factor(string, levels=chars))

# Convert a file of sites into a matrix
sites2matrix <- function(file="Hin2206sites.txt"){
    sites    <- scan(file, what="character") # read in a file of text strings
    sites    <- tolower(sites)               # make everything lowercase
    nsites   <- length(sites)                # how many sites?
    nbases   <- nchar(sites[1])              # how long is the first site?
    if(all(nchar(sites) == nbases) == T){    # return matrix, if sites all same length
        sitestr <- stringsplit(sites)        # vectorize
        goodchars <- which(sitestr %in% chars[-length(chars)]) # index ACGT pos, 
        sitestr[-goodchars] <- chars[length(chars)]            # Nify non-ACGT
        sitemat <- matrix(sitestr, nrow=nsites, ncol=nbases, byrow=TRUE) # matrix
        return(sitemat)   # output results
        }else{return(NA)} # return NA, if sites not all same length
    }

########################
# Matrix manupulations #
########################

# Position Count Matrix.  NOTE THIS HAS FIVE ROWS TO CATCH Ns
matrix2pcm   <- function(sitemat) apply(sitemat, 2, tablepos)

# Position Frequency Matrix.  NOTE THIS REMOVES THE FIFTH ROW (Ns)
pcm2pfm      <- function(pcm, pseudo=0) {
    apply(pcm[1:4,] + pseudo, 2, function(x) x/(sum(x)))
    } 
    # Note: pseudo argument is to add "pseudo-counts" if needed
    # I.e. if any cells in PCM are 0, set pseudo to 1, for example.
    # This is primarily because log(0) is undefined.
    # Other reasons in literature seem spurious, so don't bother if pcm has no zeros.

# Position-specific scoring matrix (equiprobable background frequency)
pfm2pssm <- function(observed, expected){
	expected <- matrix(0.25, dim(observed)[1], dim(observed)[2]) # weak #creates a matrix of 0.25 of the same dimentions as before
	observed * log(observed / expected, log2(length(chars)-1)) #wikipedia pssm formula
	}
	# This is the main place in these functions that could be considerably expanded.
	# Many possible scoring matrices might be used, including the pfm in the last step
	# This choice follows the Wikipedia for PSSM, as it is more flexible for non-equiprobable bases
	# An adjustment would allow accounting for GC-content
	# A further refinement would be to use the scaled values as in Mell 2012 Nar

# Make a weblogo-style matrix (equiprobable background frequency)
pfm2wblg <- function(observed, expected){
     expected <- matrix(0.25, dim(observed)[1], dim(observed)[2]) # weak
	 pssm     <- pfm2pssm(observed, expected)
	 ic       <- colSums(pssm)

   t(t(observed) * ic) # Gotta transpose!
	 }


##################
# Genome massage #
##################

# Ensure all non-ACGT are Ns for a contig
ncontigs <- function(contig, valid=chars){
    contig[-which(contig %in% valid[-length(valid)])] <- valid[length(valid)]  #valid[-length(valid)] takes first 4 chars
    return(contig)
    }

# Ensure all non-ACGT are Ns for a contig but with gaps
gapcontigs <- function(contig, valid=chars.gap){
  contig[-which(contig %in% valid[-5])] <- valid[5]  #valid[-length(valid)] takes first 4 chars
  return(contig)
}

# Factorize contig according to valid chars
fcontigs <- function(contig, valid=chars){ 
    contig <- factor(contig, levels=valid)
    return(contig)
    }


# Factorize contig according to valid chars
gapfcontigs <- function(contig, valid=chars.gap){ 
  contig <- factor(contig, levels=valid)
  return(contig)
}

# Wrapper for reading in a genome as a list
read.genome <- function(fasta){
    genome <- read.fasta(fasta) # read fasta with seqinr function
    genome <- lapply(genome, ncontigs) # non-ACGT->N
    genome <- lapply(genome, fcontigs) # reorder bases to ACGTN
    return(genome)
}

# Wrapper for reading in a genome as a list
read.genome.gaps <- function(fasta){
  genome <- read.fasta(fasta) # read fasta with seqinr function
  genome <- lapply(genome, gapcontigs) # non-ACGT->N
  genome <- lapply(genome, gapfcontigs) # reorder bases to ACGTN
  return(genome)
}

# Concatenate beginning of contig to the end to score full circle
circleit <- function(pssm, gennum) c(gennum, gennum[1:dim(pssm)[2]])

#############################
# Scoring and normalization #
#############################

# Score with a pssm using a query of the same length
scoreit <- function(pssm, query) sum(diag(pssm[query, ])) # faster way?

# Normalize a vector of scores by the best and worst sequences of a pssm
normit <- function(scores, pssm){
    best     <- apply(pssm, MARGIN=2, which.max) # best site
    worst    <- apply(pssm, MARGIN=2, which.min) # worst site
    bestest  <- scoreit(pssm, best) # score of best
    worstest <- scoreit(pssm, worst) # score of worst
    
    # Normalize scores to between 0 and 1 by worst and best sequences
    sapply(scores, function(score) (score - worstest) / (bestest - worstest))
    }

################################################
# Score a contig on both strands in all frames #
################################################

scoreContig <- function(query, scoremat, circle=T){
    
    # Motif handling
    sitelen <- dim(scoremat)[2] # length of site
    n       <- rep(0, dim(scoremat)[2]) # row for Ns
    pssmN   <- rbind(scoremat, n) # pssm with Ns
    pssmRC  <- pssmN[c((length(chars)-1):1,length(chars)), sitelen:1] #invert pssm matrix

    # Genome handing
    gennum  <- as.numeric(query) # convert contig to numeric vector
    genlen  <- length(gennum)    # length of genom
    
    # Control concatenation at end of sequence by circle parameter
    if(circle == T) gennum <- circleit(scoremat, gennum) # circularize?
    if(circle != T) gennum <- c(gennum, rep(length(chars), sitelen - 1))   # no circle?
    
    # Preallocate every site to score (time sink/memory hog?)
    genmat <- sapply(1:genlen, function(pos) gennum[pos:(pos + sitelen - 1)])

    # Score every frame, returning a matrix with scores for each strand (INEFFICIENT!!!)
    z <- apply(genmat, 2, function(subseq){
        f <- scoreit(pssmN,  subseq) # score for top strand

        r <- scoreit(pssmRC, subseq) # score for bottom strand
        return(c(f,r))
        }
        )
        # THIS IS THE TIME-KILLER.  SHOULD BE FASTER!

    z <- data.frame(t(z))      # Need to transpose for sanity
    colnames(z) <- c("w", "c") # refactor from 12 to wc
    return(z)                  # return the data.frame
    }
# Further modifications to increase speed highly desirable

################################
# Index strand with best score #
################################

# takes output from scoreContig

max.strand <- function(scorepairs){
    z <- factor(apply(scorepairs, 1, which.max)) # factorize
    levels(z) <- list(w=1, c=2) # re-level to wc instead of 12
    return(z)
    }

##########################################################################
###################### MORE FUNCTIONS GO BELOW !!! #######################
##########################################################################

# function to calculate the distance to closest USS or to closest peak
# this function compares two lists of genomic positions, the first should be the reference list. For instance a list of uptake peaks 
# the second argument should be the list that will be compared with the first list. For instance a list of USS positions
# results will always have the same length as the first list since this is the reference from which distance to the second list will be calculated
dist.USS <- function(position, USS.genome.c){                
  closest.USS<- min(abs(position - USS.genome.c))   
  return(closest.USS)                     # return result
}

#original peak finder function as described in http://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset 
argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

#original function of the peak finder: http://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset
# used to plot the uptake ratio per genomic positions as well as positions identified as peaks (red) and the smoothe line (black)
original.test <- function(w, span) {
  peaks <- argmax(x, y, w=w, span=span)
  
  plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
  lines(x, peaks$y.hat,  lwd=2) #$
  y.min <- min(y)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]),
                                    col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
}



#new peak finder function which add the cut argument.
# this argument removes peaks lower than a certain value "cut"
# This function was modified from: http://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset
# Read the link aboe to understand the span and w arguments
# x is usually a column in a dataframe with genomic positions
# y is usually a column in a dataframe with uptake ratios of x positions
# note: this function does not perfom well on entire genomes of large genomic regions (>10kb)
test <- function(x , y, w, span, cut) {
  peaks <- argmax(x, y, w=w, span=span)
  tem <- peaks[[1]] # add numbers of position to a vector instead of a list
  peaks[[4]]<- uptake$re.ratio[tem]  #get uptake ratio of positions with peaks
  trim<- which(peaks[[4]] >= cut) #eliminate positions with less than 0.5        
  pos<- peaks$x[trim]  
  ratio<- uptake$ratio[pos]
  re.ratio<- uptake$re.ratio[pos]
  peak.data<- data.frame( pos = pos, ratio = ratio, re.ratio = re.ratio)  
}  



#This function is used to plot the uptake ratio per genomic position (in grey) plus the smoothed line calculated as part of the peak finder function
# Then the function add the points identified as peaks in red. 
# This function requires to input  four parameters
# This function requires to pre-calculate p.hat (smoothed points) with the argmax function above
#IMPORTANT: This function requires the result from the function test (argument h) and the list of smoothed p.hay uptake ratio points from the function argmax see: final_peak_finder.R
peak.plot<- function(h, x, y, span = spa) {
  plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
  lines(x, p.hat,  lwd=2) #$
  y.min <- min(y)
  sapply(h$pos, function(i) lines(c(x[i],x[i]), c(y.min, p.hat[i]),
                                  col="Red", lty=2))
  points(x[h$pos], p.hat[h$pos], col="Red", pch=19, cex=1.25)
  
}



#Given a genomic position, this function retrieve the sequence +- the flanking size chosen 
# Motif retrieval function
retrieve.seq <- function(site,  genome=np, flank){
  genome[(site-flank):(site +flank)]
}


# This function creates a list including for a certain position (site): genomic site, recalculated uptake ratio
# sequence +- flanking size, maximum uptake USS score of the sequence range chosen (for each strand), how far from the chosen site is the max USS score position.  
pos.detail.2 <- function(site,  genome=np, flank){
  seq<- retrieve.seq(site = site, genome = genome, flank = flank)
  score.w<- np.uptake.scores$w[(site-flank):(site +flank)]  
  score.c<- np.uptake.scores$c[(site-flank):(site +flank)]
  max.score.w<- max(score.w)   #maximum uptake USS score of the sequence range chosen (forward strand) 
  max.score.c<- max(score.c)   #maximum uptake USS score of the sequence range chosen (reverse strand)
  max.pos.w<- (which.max(score.w) - (flank +1)) #how far from the chosen site is the max USS score position (forward strand)
  max.pos.c<- (which.max(score.c) - (flank+1))  #how far from the chosen site is the max USS score position (reverse strand)
  seq1<- paste(seq, collapse = "") # merge elements of a vector to make a one line sequence   
  re.ratio<- uptake$re.ratio[site] #recalculated uptake ratio for chosen site
  a<- list(site = site, re.ratio = re.ratio, sequence =  seq1,  max.score.w = max.score.w, max.score.c = max.score.c, position.of.max.score.w = max.pos.w, position.of.max.score.c = max.pos.c)
  
}




