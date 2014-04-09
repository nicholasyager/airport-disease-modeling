# SIR_Means.R
# Nicholas A. Yager and Matthew Taylor.
#
# This script plots the mean SIR response for the disease propagating through the network
#

filelist = dir("sir", full.names=T)

for (file_num in 1:length(filelist)) {
  file <- read.csv(filelist[file_num])
  
  S <- file$s
  I <- file$i
  R <- file$r
  
  for ( i in 1:(40-length(S))) {
    S = append(S,S[length(S)])
  }
  
  for ( i in 1:(40-length(R))) {
    R = append(R,R[length(R)])
  }
  
  for ( i in 1:(40-length(I))) {
    I = append(I,I[length(I)])
  }
  
  if (file_num != 1) {
    
    allS <- cbind(allS, S)
    allI <- cbind(allI,I)
    allR <- cbind(allR,R)
    
  } else {
    
    allS <- data.frame(S)
    allI <- data.frame(I)
    allR <- data.frame(R)
    
  }
  
}

meanS <- apply(as.matrix(allS), 1, mean)
meanI <- apply(as.matrix(allI), 1, mean)
meanR <- apply(as.matrix(allR), 1, mean)

stddevS <-apply(as.matrix(allS), 1, sd)
stddevI <-apply(as.matrix(allI), 1, sd)
stddevR <-apply(as.matrix(allR), 1, sd)

totals = meanS + meanI + meanR

plot(0,type="n", ylim=c(0,3300), xlim=c(0,30))

apply(as.matrix(allS), 2, function(x) points(x,col="blue",pch=1,cex=0.5,))
apply(as.matrix(allI), 2, function(x) points(x,col="green",pch=1,cex=0.5))
apply(as.matrix(allR), 2, function(x) points(x,col="red",pch=1,cex=0.5))

plot(meanS/totals,
     main="Mean SIR Response of Infection in Weighted Airport Network",
     xlab="Time (ticks)", ylab="Proportion of Population",
     type="l",col="blue", ylim=c(0,1),  lwd=2)
lines(meanI/ totals, col="green", lty="dashed", lwd=2)
lines(meanR/ totals, col="red", lty="dotted", lwd=2)

