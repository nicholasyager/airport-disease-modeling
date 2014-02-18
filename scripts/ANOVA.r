library(plotrix)

R <- as.data.frame(read.csv("random.matrix",header=F))
D <- as.data.frame(read.csv("degree.matrix",header=F))
B <- as.data.frame(read.csv("betweenness.matrix",header=F))

r <- data.frame()
k = 1
for (i in c(1,2,3,4,5)) {
  for (j in 1:50) {
    r[k,1] <- "random"
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- R[j,i]
    k = k + 1
  }
}
for (i in c(1,2,3,4,5)) {
  for (j in 1:50) {
    r[k,1] <- "degree"
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- D[j,i]
    k = k + 1
  }
}
for (i in c(1,2,3,4,5)) {
  for (j in 1:50) {
    r[k,1] <- "betweenness"
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- B[j,i]
    k = k + 1
  }
}

data.stacked <- r
names(data.stacked) <- c("qtype","effort","infected")

#data.stacked = stack(data)

data.aov = aov(data.stacked$infected~data.stacked$qtype*data.stacked$effort)
summary(data.aov)
TukeyHSD(data.aov)

M <- tapply(data.stacked$infected,list(data.stacked$qtype, data.stacked$effort),median)
M <- M[,order(M[2,],decreasing=T)]

CI95 = qt(1-0.05/2,9)*(sd(data.stacked$infected)/sqrt(length(data.stacked$infected)))

bp = barplot(M,beside=T, 
             xlab="Quarantine Effort (% of airports closed)", ylab="Number of infected airports",
             main="Comparison of quarantine strategies",
             ylim=c(0,4500))
legend(x=16.1,y=4500, legend=c("Betweenness","Degree","Random"), fill=c("gray30","gray60","gray90"), cex=0.75)
plotCI(bp, M, CI95, add=T, pch=NA)
text(x=as.vector(bp),y=as.vector(M)+CI95,c("A","A","A",
                                           "A","A","A",
                                           "B","C","A",
                                           "D","D","A",
                                           "D","D","A"), pos=3)

abline(h=0)
