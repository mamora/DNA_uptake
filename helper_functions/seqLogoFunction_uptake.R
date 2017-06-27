
# WISH LIST: 
#	1. Only plot the consensus base, without the proportional stack.
#	2. Possibly simpler: Let the -ic bases go negative!

pwm2ic2 <- function (pwm) 
{
    npos <- ncol(pwm)
    base <- integer(length = npos)	# this handles excluded/absent bases
    for (i in 1:npos) {
    	base[i] <- length(which(pwm[, i] != 0))
    	}
    ic <- numeric(length = npos)
    for (i in 1:npos) {
        ic[i] <- log2(base[i]) + sum(sapply(pwm[, i], function(x) {
            if (x > 0) {
                x * log2(x)
            } else {
                0
            }
        }))
    }
    ic
}



seqLogo2 <- function (pwm, ic.scale = TRUE, 
                      xaxis = TRUE, yaxis = TRUE, 
                      xfontsize = 12, yfontsize = 12, 
                      maradj = 0, fontadj = 3, yscale = 2
                      ) 
{
    if (class(pwm) == "pwm") {
        pwm <- pwm@pwm
    }
    else if (class(pwm) == "data.frame") {
        pwm <- as.matrix(pwm)
    }
    else if (class(pwm) != "matrix") {
        stop("pwm must be of class matrix or data.frame")
    }
    if (any(abs(1 - apply(pwm, 2, sum)) > 0.01)) 
        stop("Columns of PWM must add up to 1.0")
    chars <- c("A", "C", "G", "T")
    letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos <- ncol(pwm)
    if (ic.scale) {
        ylim <- yscale
        ylab <- "bits"
        facs <- pwm2ic2(pwm)
    }
    else {
        ylim <- 1
        ylab <- "probability"
        facs <- rep(1, npos)
    }
    wt <- 1
    x.pos <- 0
    for (j in 1:npos) {
        column <- pwm[, j]
        hts <- 0.95 * column * facs[j]
        letterOrder <- order(hts)
        y.pos <- 0
        for (i in 1:4) {
            letter <- chars[letterOrder[i]]
            ht <- hts[letterOrder[i]]
            if (ht > 0) 
                letters <- seqLogo:::addLetter(letters, letter, x.pos, 
                  y.pos, ht, wt)
            y.pos <- y.pos + ht + 0.01
        }
        x.pos <- x.pos + wt
    }
    grid.newpage()
    bottomMargin = ifelse(xaxis, maradj + xfontsize/3.5, maradj)
    leftMargin = ifelse(yaxis, maradj + yfontsize/3.5, maradj)
    pushViewport(plotViewport(c(bottomMargin, leftMargin, maradj, 
        maradj)))
    pushViewport(dataViewport(0:ncol(pwm), 0:ylim, name = "vp1"))
    grid.polygon(x = unit(letters$x, "native"), y = unit(letters$y, 
        "native"), id = letters$id, gp = gpar(fill = letters$fill, 
        col = "transparent"))
    if (xaxis) {
        grid.xaxis(at = seq(0.5, ncol(pwm) - 0.5), label = 1:ncol(pwm), 
            gp = gpar(fontsize = xfontsize-fontadj))
        grid.text("position", y = unit(-3, "lines"), gp = gpar(fontsize = xfontsize))
    }
    if (yaxis) {
        grid.yaxis(gp = gpar(fontsize = yfontsize-fontadj))
        grid.text(ylab, x = unit(-3, "lines"), rot = 90, gp = gpar(fontsize = yfontsize))
    }
    popViewport()
    popViewport()
    par(ask = FALSE)
}

#ls()
#save(pcm, pinrn, pwm2ic2, seqLogo2, file="Outputs/functionsMarch11.R")