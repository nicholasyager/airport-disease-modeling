M <- matrix(nrow=330, ncol=330, byrow=T)

for (j in 1:330) {
  
  for (i in 1:330) {
    sus = 330 - i
    # Distance from 0,o
    d = sqrt((j)^2+(i)^2)
    M[j,i] <-  -(i+j)^2
    
  }
  
}

image(M, col=rainbow(100))

filelist = dir()

for (generation in 1:length(filelist)-1){
  print(generation)
 
  pareto = read.csv(paste("pareto_",generation,".csv",sep=""))
  
  colFunc <- colorRampPalette(c("lightgray", "black"))
  colFunc2 <- colorRampPalette(c("blue", "red"))
  colors = colFunc(50)
  colors2 = rainbow(100)

#   image(M, col=colors2,
#         main=paste("Fitness of Quarantine Strategies - Generation",generation),
#         xlab="Proportion Airports Quarantined",
#         ylab="Proportion Airports Infected")
  par(mfrow=c(1,2))
  image(M, col=rainbow(100),
         main=paste("Fitness of Quarantine Strategies - Generation",generation),
         xlab="Airports Quarantined",
         ylab="Airports Infected")
  points(pareto$closures/3307, pareto$infected/3307,
        col=colors ,
        pch=16)
       
  hist(pareto$fitness)
  
  abline(v=mean(pareto$fitness), col="red")

}