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

types = c("betweenness","random","cluster")

par(mfrow=c(2,2))
for (type_num in 1:length(types)){
  
  type = types[type_num]
  
  random <- read_to_M(dir(type,full.names=T))
  medians = apply(random, 2, mean)
  
  effort = seq(0, 100, 5)
  data = as.data.frame(cbind(effort, medians))
  names(data) <- c("Effort", "Median")
  write.csv(data, paste(type,".csv",sep=""), row.names=F)
  write.matrix(random, paste(type,".matrix",sep=""), sep=",")
  
  # A nice distribution graph -------------------------------------------
  
  boxplot(random,
          main=paste(type,"-based quaratine strategy",sep=""),
          xlab="Quarantine Effort", ylab="Proportion of airports infected",
          xaxt = "n")
  
  axis(1,1:21, seq(0,100,5))
  lines(apply(random, 2, mean), col="red")
  lines(apply(random, 2, median), col="blue")
  
}
par(mfrow=c(1,1))
#image(random, col=rainbow(3000))
# 
# 
#  barplot(medians/3307,
#          names=effort,
#          main="Effect of quarantine strategies on median number of infections.",
#          xlab="Degree-based quarantine effort",
#          ylab="Proportion of airports infected",
#          ylim=c(0,1))
#  abline(h=0)
