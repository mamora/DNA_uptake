# logo work
#source("http://bioconductor.org/biocLite.R")
#biocLite("seqLogo")


setwd("C:/Users/marcelo/Dropbox/ Uptake ms/DegenerateStuff/")

library(seqLogo)
library(RColorBrewer)

#### Contains a bunch of functions (probably hacky)
# Main one is SeqLogo2, which is my hack to get the right object 
load("Outputs/functionsMarch11.R")

inpPCM <- read.csv("Outputs/inppcm.csv", row.names=1)
perPCM <- read.csv("Outputs/perpcm.csv", row.names=1)

# This gets the consensus of the right 31 bases and makes a numeric string too.
consensus <- "ATGCCAAAGTGCGGTTAATTTTTACAGTATTTTTGGGTTCGA"
consensus <- substring(consensus, 6, 36)
connum    <- as.numeric(factor(unlist(strsplit(consensus, ""))))
degen     <- 0.24 # This is the predicted dengeneracy of the input per position
sampn     <- 1e7  # This is the subsample size to make both datasets the same size.

# inpcons <- "AAAGTGCGGTTAATTTTTACAGTATTTTTGG"
# connum
# as.numeric(factor(unlist(strsplit(inpcons, ""))))

# subset looks awesome
# so does extract

# Expected Input... hacky as can be.  2010.
expPCM           <- matrix(degen/3, length(connum), 4)
rownames(expPCM) <- rownames(inpPCM)
colnames(expPCM) <- colnames(inpPCM)
for(i in 1:31){
	expPCM[i, connum[i]] <- 0.76
	}
expPCM <- sampn * expPCM

# WebLogos

# This is how the input differs from expected
z <- inpPCM/expPCM / rowSums(inpPCM/expPCM)
seqLogo2(t(z))
sum(pwm2ic2(t(z))) # bits inp vs exp


# This is the final derived motif,
# I.e. this is the "uptake motif" based on the degenerate uptake experiment.
# Essentially, this says what the motif would look like IF the bases had been equiprobable at each position,
# It accounts for the position-specific degeneracy.
y <- perPCM/inpPCM / rowSums(perPCM/inpPCM)
seqLogo2(t(y))
sum(pwm2ic2(t(y))) # bits per vs inp

# Thisis how the motif would look, if we'd just assumed the input degeneracy was exact
x <- perPCM/expPCM / rowSums(perPCM/expPCM)
seqLogo2(t(x))

#########################################
###### STOP #############################
# Below here is probably not useful #

# ########################
# # Matches to consensus #
# ########################
# 
# ##### New
# 
# pertab <- read.delim("~/Desktop/degenerate/hist/perip.miss.count.hist", header=F)
# inptab <- read.delim("~/Desktop/degenerate/hist/input.miss.count.hist", header=F)
# 
# sum(pertab$V2)
# sum(inptab$V2)
# 
# barplot(rbind(pertab$V3, inptab$V3[1:length(inptab$V3)-1]), names.arg=0:18, beside=T)
# 
# 
# barplot(rbind(pertab$V3, inptab$V3[1:19], (dbinom(0:31, 31, 0.24)*1e7)[1:19]), names.arg=0:18, beside=T)
# 
# 
# dbinom(0, 31, 0.24) * 1e7
# 
# # mismatch histogram
# quartz(, 4, 4, file="../mismatchDist.png", type="png", dpi=300)
# plot(0:18, (dbinom(0:31, 31, 0.24)*1e7)[1:19], type="o", col="grey", ylim=c(0, 2e6),
#      ylab="millions of fragments", xlab="# of non-consensus bases in fragment",
#      axes=F
#      )
# points(0:18, inptab$V3[1:19], type="o", col="blue")
# points(0:18, pertab$V3[1:19], type="o", col="red")
# axis(1, at=0:18, cex=0.8, las=2)
# axis(2, at=seq(0, 2e6, by=5e5), labels=seq(0, 2, by=0.5), las=1)
# legend("topright", legend=c("expected", "input", "recovered"), cex=0.8, col=c("grey", "blue", "red"), pch=21, bty="n")
# dev.off()
# 
# quartz(, 4, 4, file="../mismatchFrac.png", type="png", dpi=300)
# barplot(pertab$V3*0.079/inptab$V3[1:19], names.arg=0:18, ylim=c(0,1), cex.names=0.8, las=2,
#         ylab="fraction of input taken up", xlab="# of non-consensus bases in fragment"
#         )
# #abline(h=0.079, col="red")
# dev.off()
# 
# #############
# # Diversity #
# #############
# 
# perdiv <- read.delim("~/Desktop/degenerate/hist/perip.uniq.count.hist", header=F)
# inpdiv <- read.delim("~/Desktop/degenerate/hist/input.uniq.count.hist", header=F)
# 
# 
# perdiv <- perdiv[dim(perdiv)[1]:1, ]
# inpdiv <- inpdiv[dim(inpdiv)[1]:1, ]
# 
# 
# # "histograms"
# #quartz(,4,4)
# quartz(,4,4, file="../DiversInput.png", type="png", dpi=300)
# par(mar=c(4,4,1,1))
# plot(inpdiv, type="h", log="xy", col="blue",
#      xlim=c(1, 1e6), xlab="# of occurrences of a variant",
#      ylim=c(1, 1e7), ylab="# of variants"
#      )
# dev.off()
# #quartz(,4,4)
# quartz(,4,4, file="../DiversPerip.png", type="png", dpi=300)
# par(mar=c(4,4,1,1))
# plot(perdiv, type="h", log="xy", col="red",
#      xlim=c(1, 1e6), xlab="# of occurrences of a particular variant",
#      ylim=c(1, 1e7), ylab="# of variants"
#      )
# dev.off()
# 
# # Another way
# quartz(,4,4)
# par(mar=c(4,4,1,1))
# plot(perdiv, pch=21, cex=0.5, lwd=0.5, log="xy", col="salmon",
#      xlim=c(1, 1e6), xlab="# of occurrences of a particular variant",
#      ylim=c(1, 1e7), ylab="# of variants", axes=F, cex.lab=0.8
#      )
# points(inpdiv, pch=22, cex=0.5, col="lightblue", lwd=0.8)
# axis(1, at=10^(0:6), labels=10^(0:6), las=2, cex.axis=0.8)
# axis(2, at=10^(0:7), labels=10^(0:7), las=2, cex.axis=0.8)
# 
# # histograms of histograms?
# 
# #percut <- cut(perdiv$V1, breaks=2^seq(0, 20), right=F)
# #inpcut <- cut(inpdiv$V1, breaks=2^seq(0, 20), right=F)
# 
# percut <- cut(perdiv$V1, breaks=10^seq(-1, 6), right=T)
# inpcut <- cut(inpdiv$V1, breaks=10^seq(-1, 6), right=T)
# 
# 
# perbin <- tapply(perdiv$V2, INDEX=percut, FUN=sum)
# inpbin <- tapply(inpdiv$V2, INDEX=inpcut, FUN=sum)
# 
# quartz(, 4, 4)
# quartz(,4,4, file="../DiversBinned.png", type="png", dpi=300)
# 
# par(mar=c(5, 5, 1, 1))
# z <- barplot(rbind(inpbin, perbin), log="y", las=2, beside=T, axes=F, names.arg=rep("", length(perbin)), 
#              xlab="binned occurrences (upper-bound)",
#              ylab="# of sequences", cex.lab=0.8, cex.axis=0.8, ylim=c(0.1, 1e7),
#              col=c("lightblue", "salmon"), border=F
#              )
# axis(1, at=(z[1,]+z[2,])/2, labels=10^(0:6), las=2, cex.axis=0.8)
# #axis(1, at=(z[1,]+z[2,])/2, labels=2^(0:19), las=2, cex.axis=0.8)
# #axis(1, at=(z[1,]+z[2,])/2, labels=2^seq(1:20)-1, las=2, cex=0.8)
# axis(2, at=c(0.1, 10^(0:7)), labels=c(0, 10^(0:7)), las=1, cex.axis=0.8)
# legend("topright", legend=c("input", "recoverd"), col=c("lightblue", "salmon"), pch=15, cex=0.8, bty="n")
# dev.off()
# 
# #z <- barplot(log2(perbin/inpbin), las=2, axes=F, names.arg=F)
# #axis(1, at=z, labels=2^seq(1:20)-1, las=2, cex=0.8)
# #axis(2, at=-1:5, labels=2^(-1:5), las=1)
# #abline(h=0)
# 
# 
# 
# 
# 
# 
# 
# ##### WHY IS THIS DIFFERENT THAN PREVIOUS????? ######
# 
# 
# ?factorial
# 
# # possible classes of mismatch
# n <- 31
# k <- 0
# 3^k * choose(n, k)
# 
# factorial(31)/(factorial(29)*factorial(2)) == choose(n, k)
# # why false?
# 
# binomExp <- vector()
# for(i in 1:19){
# 	binomExp[i] <- 3^(i-1) * choose(n, (i-1))
# 	}
# binomExp
# 
# 
# 
# 
# # Answer:  didn't trim 1st degenerate base
# # didn't help!  because this was trimmed in the recent iteration...  problem with Levenshtein edit distance?
# # main problem is that there are too many classes of 2-off in both matrices
# # equation? nope binomial coefficient times number of ways to have x mismatches...  should be right.
# 
# #### Oldie
# 
# perdist <- scan("~/Dropbox/Uptake/Outputs/perdist.txt")
# inpdist <- scan("~/Dropbox/Uptake/Outputs/inpdist.txt")
# 
# length(perdist)
# length(inpdist)
# 
# rev(tabulate(perdist+1, 32))
# rev(tabulate(inpdist+1, 32))
# 
# 
# 
# 
# barplot(dbinom(0:31, 31, 0.76), names.arg=0:31)
# 
# cbind(0:32, dbinom(0:31, 31, 0.76))
# 
# plot(0:31, dbinom(0:31, 31, 0.76)*1e7, type="l", lwd=2, col="grey", xlab="# of matches to consensus USS", ylab="Frequency of Reads")
# points(0:31, rev(tabulate(inpdist+1, 32)), type="l", lwd=2, col="red")
# points(0:31, rev(tabulate(perdist+1, 32)), type="l", lwd=2, col="blue")
# 
# 
# 
# which.max(dbinom(0:31, 31, 0.76)*1e7)
# which.max(rev(tabulate(inpdist+1, 32)))
# which.max(rev(tabulate(perdist+1, 32)))
# 
# 
# expect <- dbinom(0:31, 31, 0.76)*1e7
# input  <- rev(tabulate(inpdist+1, 32))
# perip  <- rev(tabulate(perdist+1, 32))
# 
# plot(0:31, expect, type="l", lwd=1, col="grey",
#      main="Distribution of matches to consensus", 
#      xlab="# of matches to 31 bp consensus", ylab="millions of reads", 
#      ylim=c(0, 2e6), xlim=c(15, 31), xaxt="n", yaxt="n"
#      )
# axis(1, at=0:31, cex.axis=0.8)
# axis(2, at=seq(0, 2e6, by=5e5), labels=c(0.0, 0.5, 1.0, 1.5, 2.0), las=1)
# 
# #points(0:31, expect, type="h", lwd=2, col="grey")
# points(0:31, expect, pch=20, col="grey")
# 
# points(0:31, input, type="l", lwd=1, col="red")
# #points(0:31, input, type="h", lwd=2, col="red")
# points(0:31, input, pch=20, col="red")
# 
# points(0:31, perip, type="l", lwd=1, col="blue")
# #points(0:31, perip, type="h", lwd=2, col="blue")
# points(0:31, perip, pch=20, col="blue")
# 
# legend("topleft", legend=c("expected", "input", "recovered"), pch=20, col=c("grey", "red", "blue"), bty="n")
# 
# 
# 
# barplot(dbinom(0:31, 31, 0.76)*1e7)
# , rev(tabulate(inpdist+1, 32)))
# 
# 
# barplot(rev(tabulate(inpdist+1, 32)), add=T)
# , rev(tabulate(inpdist+1, 32)), rev(tabulate(perdist+1, 32)), beside=T, names.arg=0:31)
# 
# 
# 
# ?plot