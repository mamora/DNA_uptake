

#original peakfinder script
argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}


#modified peakfinder script
peakfinder <- function(x , y, w, span, cut) {
  peaks <- argmax(x, y, w=w, span=span)
  tem <- peaks[[1]] # add numbers of position to a vector instead of a list
  peaks[[4]]<- uptake$re.ratio[tem]  #get uptake ratio of positions with peaks
  trim<- which(peaks[[4]] >= cut) #eliminate positions with less than 0.5        
  pos<- peaks$x[trim]  
  ratio<- uptake$re.ratio[pos]
  peak.data<- data.frame( pos = pos, ratio = ratio)  
}

# split uptake dataframe and extract positions
split_pos<- function(uptake){
overlap<- 1999
diff<- (length(uptake$pos) - 10000*191) #calculating the length of the last chunk 
last.chunck<- (10000 + (length(uptake$pos) - 10000*191)) #making a vector of the size of the last chunck
p<- rep(10000,190) #make a vector of 19 times 100000
p<- c(p,(last.chunck))  # make a vector of the all the chunck sizes
pos.l<- split(uptake$pos, rep(1:191, p)) #split the positions of the genome in one list
pos.l2<- pos.l

#make a list of the split genome with 2000bp overlaps                
for (i in 2:190){
  a<- i + 1
  b<- i - 1
  pos.l2[[i]]<- c((pos.l[[b]][10000]- overlap):(pos.l[[b]][10000]),pos.l[[i]],(pos.l[[a]][1]:(pos.l[[a]][1]+ overlap)))
}
#circularize
pos.l2[[1]]<- c((pos.l[[191]][14490]- overlap):(pos.l[[191]][14490]),pos.l[[1]],(pos.l[[2]][1]:(pos.l[[2]][1]+ overlap)))
pos.l2[[191]]<- c((pos.l[[190]][10000]- overlap):(pos.l[[190]][10000]),pos.l[[191]],(pos.l[[1]][1]:(pos.l[[1]][1]+ overlap)))
#end of function
return(pos.l2)
}

# split uptake dataframe and extract ratios
split_ratio<- function(uptake){
  overlap<- 1999
  diff<- (length(uptake$pos) - 10000*191) #calculating the length of the last chunk 
  last.chunck<- (10000 + (length(uptake$pos) - 10000*191)) #making a vector of the size of the last chunck
  p<- rep(10000,190) #make a vector of 19 times 100000
  p<- c(p,(last.chunck))  # make a vector of the all the chunck sizes
  ratio.l<- split(uptake$re.ratio, rep(1:191, p)) # split the ratios of the genome in another list
  ratio.l2<- ratio.l

  for (i in 2:190){
  a<- i + 1
  b<- i - 1
  ratio.l2[[i]]<- c(ratio.l[[b]][(10000 - overlap):10000],ratio.l[[i]],(ratio.l[[a]][1:(1+ overlap)]))
}

#circularize
ratio.l2[[1]]<- c(ratio.l[[191]][(14490 - overlap):14490],ratio.l[[1]],(ratio.l[[2]][1:(1+ overlap)]))
ratio.l2[[191]]<- c(ratio.l[[190]][(10000 - overlap):10000],ratio.l[[191]],(ratio.l[[1]][1:(1+ overlap)]))

return(ratio.l2)
#end of function
}


               





