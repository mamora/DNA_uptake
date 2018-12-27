#########################################################################################
#         Functions to make plots of uptake maps of observed vs predicted model         #
#########################################################################################

x.axis.num<- seq(start.pos, end.pos, by = spacing_xaxis) 

################################
# add arrows as USSs function  #
################################
# function to add arrows in places where there are USS in a plot
make.arrows.plot<- function(p0, start, endf){
  data1<- USS.scores[(start):(endf),]
  f<- which(data1$max >= arrow.uss.cut)
  t2 <- data1[f,]
  st<- c()
  end<- c()
  for (j in 1:length(t2$w)){
    if (t2$strand[j] == "w"){
      st[j]<- t2$V1[j]
      end[j]<- t2$V1[j] + 30
    }else{
      st[j]<- t2$V1[j] + 29
      end[j]<- t2$V1[j] - 1
    }
  }
  segment_data = data.frame(
    x = st,
    xend = end, 
    y = arrows.y.height, # place where I want arrows to be in the y axis
    yend = arrows.y.height 
  )
  p0<- p0 + geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),
                         arrow = arrow(length = unit(0.2, "cm"), type = "closed"), colour = "red")
  return(p0)
}

make.arrows2.plot<- function(p0, start, endf){
  data1<- USS.scores[(start):(endf),]
  f<- which(data1$max >= arrow_2.uss.cut & data1$max < arrow.uss.cut  )
  t2 <- data1[f,]
  st<- c()
  end<- c()
  for (j in 1:length(t2$w)){
    if (t2$strand[j] == "w"){
      st[j]<- t2$V1[j]
      end[j]<- t2$V1[j] + 30
    }else{
      st[j]<- t2$V1[j] + 29
      end[j]<- t2$V1[j] - 1
    }
  }
  segment_data = data.frame(
    x = st,
    xend = end, 
    y = arrows.y.height, # place where I want arrows to be in the y axis
    yend = arrows.y.height 
  )
  p0<- p0 + geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),
                         arrow = arrow(length = unit(0.2, "cm"), type = "closed"), colour = "blue")
  return(p0)
}


####################################################
# collection of different themes to use in a plot  #
####################################################

theme_background_grey<- theme(plot.margin = unit(c(1,1,1,1),"cm"),
                              legend.position = "bottom",
                              panel.grid.minor = element_line(colour="white", size=0.5),
                              plot.title = element_text(size = title.size, face = "bold", hjust = 0.5),
                              axis.text  = element_text(size= xy.num.text),
                              axis.title = element_text(size = xy.title.size, face = "bold")) 


theme_background_white <- theme_bw() + 
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = title.size, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= xy.num.text),
        axis.title = element_text(size = xy.title.size, face = "bold"))

theme_background_minimal <- theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1),"cm"),
        legend.position = "bottom",
        panel.grid.minor = element_line(colour="white", size=0.5),
        plot.title = element_text(size = title.size, face = "bold", hjust = 0.5),
        axis.text  = element_text(size= xy.num.text),
        axis.title = element_text(size = xy.title.size, face = "bold"))

theme_background_translucid <- theme(panel.border = element_blank(),
                                     legend.key = element_blank(),
                                     axis.ticks = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.text.x = element_blank(),
                                     panel.grid = element_blank(),
                                     panel.grid.minor = element_blank(), 
                                     panel.grid.major = element_blank(),
                                     panel.background = element_blank(),
                                     plot.background = element_rect(fill = "transparent",
                                                                    colour = NA))

################################################################
# Function to plot uptake map with SD of replicates            #
################################################################

standard.sd.plot <- function(data = simple_model, theme = theme_background_grey,  add.arrows = TRUE){
  
  p<- data %>% dplyr::filter(pos >= start.pos & pos <= end.pos) %>% 
    ggplot() +
    geom_line(aes(x = pos, y = uptake)) +
    geom_ribbon(aes(x = pos, ymin=uptake-sd_short, ymax=uptake+sd_short), color = "blue", alpha =0.3) +
    scale_x_continuous(breaks = x.axis.num, expand = c(0, 0)) +
    scale_y_continuous(limits = c(y.axis.start, y.axis.height), expand = c(0, 0))+
    guides(colour = guide_legend(override.aes = list(size=legend.size))) 
  
  temp<- p + theme
  if (add.arrows == TRUE){
    temp<- p + theme
    p1<- make.arrows.plot(p0 = temp, start =start.pos, endf = end.pos)
    return(p1)
  } else {
    temp<- p + theme
    return(temp)
  }
}

add_more_arrows<-function(plot = t){
  p1<- make.arrows2.plot(p0 = plot, start =start.pos, endf = end.pos)
  return(p1)
}  



################################################################
# Function to plot uptake map of observed vs predicted uptake  #
################################################################
standard.plot <- function(data = simple_model, theme = theme_background_grey,  add.arrows = TRUE){
  
  p<- data %>% dplyr::filter(pos >= start.pos & pos <= end.pos) %>% 
    ggplot() +
    geom_point(aes(x = pos, y = uptake), shape = 20, size = 1) +
    scale_x_continuous(breaks = x.axis.num, expand = c(0, 0)) +
    scale_y_continuous(limits = c(y.axis.start, y.axis.height), expand = c(0, 0))+
    guides(colour = guide_legend(override.aes = list(size=legend.size))) 
  
  temp<- p + theme
  if (add.arrows == TRUE){
    temp<- p + theme
    p1<- make.arrows.plot(p0 = temp, start =start.pos, endf = end.pos)
    return(p1)
  } else {
    temp<- p + theme
    return(temp)
  }
}

add_more_arrows<-function(plot = t){
  p1<- make.arrows2.plot(p0 = plot, start =start.pos, endf = end.pos)
  return(p1)
}
#############################################################################################
#            Function to plot uptake map with  scientific notation xaxis                    #
#############################################################################################
scientif.plot <- function(data = simple_model, theme = theme_background_grey,  add.arrows = TRUE){
  p<- data %>% dplyr::filter(pos >= start.pos & pos <= end.pos) %>% 
    ggplot() +
    geom_point(aes(x = pos, y = uptake), shape = 20, size = 1) +
    scale_x_continuous(breaks = x.axis.num, expand = c(0, 0), 
                       labels=function(n){format(n, scientific = TRUE, digits = 3)})+
    scale_y_continuous(limits = c(y.axis.start, y.axis.height), expand = c(0, 0))+
    labs(x = x.axis.title, y = y.axis.title) +
    ggtitle(figure.title) +
    guides(colour = guide_legend(override.aes = list(size=legend.size))) 
  
  temp<- p + theme
  if (add.arrows == TRUE){
    temp<- p + theme
    p1<- make.arrows.plot(p0 = temp, start =start.pos, endf = end.pos)
    return(p1)
  } else {
    temp<- p + theme
    return(temp)
  }
}

##########################################################################################
#                  Function to plot uptake map using log2 of the y axis                  #
##########################################################################################

log2.plot <- function(data = simple_model, theme = theme_background_grey,  add.arrows = TRUE){ 
  
  log.p<- data %>% dplyr::filter(pos >= start.pos & pos <= end.pos) %>% 
    ggplot() +
    geom_point(aes(x = pos, y = log2(uptake)), shape = 20, size = 1) +
    scale_x_continuous(breaks = x.axis.num, expand = c(0, 0), 
                       labels=function(n){format(n, scientific = TRUE, digits = 4)})+
    scale_y_continuous(limits = c(-13, 7), breaks = c(-13,-11,-9,-7,-5,-3,-1,1,3,5,7), 
                       labels = c(0.0001,0.0005,0.002,0.008,0.031,0.125,0.5,2,8,32,128), 
                       expand = c(0, 0))+
    labs(x = x.axis.title, y = y.axis.title) +
    ggtitle(figure.title) +
    guides(colour = guide_legend(override.aes = list(size=legend.size))) 
  
  
  if (add.arrows == TRUE){
    temp<- log.p + theme
    p1<- make.arrows.plot(p0 = temp, start =start.pos, endf = end.pos)
    return(p1)
  } else {
    temp<- log.p + theme
    return(temp)
  }
}
