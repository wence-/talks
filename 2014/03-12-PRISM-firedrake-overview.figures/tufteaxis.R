`tufteaxis` <- function(side, amin=NULL, amax=NULL, at=NULL, mingap=0.5, digits=2, col="black", lwd=0.7) {

  if ( side == 1 || side == 3 ) {
    islog <- par("xlog")
  } else {
    islog <- par("ylog")
  }

  if ( is.null(at) ) {
    ticks <- axTicks(side)
  } else {
    ticks <- at
  }

  domin <- !is.null(amin)
  domax <- !is.null(amax)

  if ( !domin )
    amin <- ticks[1]

  if ( !domax )
    amax <- ticks[length(ticks)]
  
  ticks <- ticks[(ticks >= amin) & (ticks <= amax)]

  nticks <- length(ticks)

  if ( islog ) {
    gap <- (log10(ticks[nticks]) - log10(ticks[nticks-1])) * mingap
  } else {
    gap <- (ticks[nticks] - ticks[nticks-1]) * mingap
  }

  
  nticks <- length(ticks)
  firsttick <- ticks[1]
  lasttick <- ticks[nticks]
  
  # If max tick will be too close to the last tick, replace it,
  #  otherwise append it
  if ( domax ) {
    if (islog && (log10(amax) - log10(lasttick) < gap)) {
      ticks[nticks]<-amax
    } else if (amax - lasttick < gap) {
      ticks[nticks]<-amax	
    } else {	
      ticks<-c(ticks,amax)
    }
  }
  
  # Similarly for first tick
  if ( domin ) {
    if (islog && (abs(log10(amin)-log10(firsttick)) < gap)) {
      ticks[1]<-amin	
    } else if (firsttick - amin < gap) {
      ticks[1]<-amin	
    } else {	
      ticks<-c(amin, ticks)
    }
  }
  
  nticks <- length(ticks)

  if ( domin )
    lmin <- format(round(ticks[1], digits), nsmall=digits, trim=TRUE)
  else
    lmin <- format(ticks[1], trim=TRUE)
  if ( domax )
    lmax <- format(round(ticks[nticks], digits), nsmall=digits, trim=TRUE)
  else
    lmax <- format(ticks[nticks], trim=TRUE)
  

  mid <- format(ticks[2:(nticks-1)], trim=TRUE)

  labels = c(lmin, mid, lmax)

  bg <- par("bg")
  if ( bg == "transparent" ) {
    bg <- "white"
  }

  if ( domin || domax ) {
    ## Draw the axis and tickmarks
    axis(side, ticks, labels, col=col, lwd=lwd)
  } else {
    ## Just use the default if no min and max are supplied
    axis(side, ticks, col=col, lwd=lwd)
  }
}
                        
