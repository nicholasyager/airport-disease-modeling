library(MASS)

filelist = dir()


read_to_M <- function(filelist) {
  print(filelist)
  toremove = c()
  for (i in 1:length(filelist)) {
    filesize  = file.info(filelist[i])$size
    print(filesize)
    if (filesize == 0) {
      toremove = append(toremove, i)
    }
  }
  if (length(toremove) > 0) {
    filelist = filelist[-toremove]
  }
  
  number_of_strategies = length(filelist)
  
  M = matrix(ncol = 21, nrow = number_of_strategies, byrow=T)
  
  # Fill the matrix with data
  
  for( strategy in 1:number_of_strategies ) {
    
    data = read.csv(filelist[strategy])
    
    for (row in 1:length(data$effort)) {
      M[strategy,row] <- data$total_infected[row]
    }
    
  }
  return(M)
}


random <- read_to_M(dir("random",full.names=T))

# Each Row is a data set for an effort. use accordingly.

# 10% effort is row 11

medians = apply(random, 2, median)

effort = seq(0, 100, 5)
plot(effort, medians/3308,
     main="Effect of quarantine strategies on median number of infections",
     xlab="Betweenness-based quarantine effort (%)",
     ylab="Proportion of airports infected", pch=19, type="b",
     ylim=c(0,1))

abline(h=0, lty="dotted")

data = as.data.frame(cbind(effort, medians/3308))
names(data) <- c("Effort", "Median")



write.csv(data, "random.csv", row.names=F)
write.matrix(random, "random.matrix", sep=",")

image(random, col=rainbow(3000))
# 
# 
#  barplot(medians/3307,
#          names=effort,
#          main="Effect of quarantine strategies on median number of infections.",
#          xlab="Degree-based quarantine effort",
#          ylab="Proportion of airports infected",
#          ylim=c(0,1))
#  abline(h=0)
