library(plotrix)
library(agricolae)

R <- as.data.frame(read.csv("random.matrix",header=F))
B <- as.data.frame(read.csv("betweenness.matrix",header=F))
C <- as.data.frame(read.csv("cluster.matrix",header=F))

r <- data.frame()
k = 1

efforts = c(1,3,7,11,17)

for (i in efforts) {
  for (j in 1:length(R[,1])) {
    r[k,1] <- "Random"
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- R[j,i]
    k = k + 1
  }
}
for (i in efforts) {
  for (j in 1:length(B[,1])) {
    r[k,1] <- "Betweenness Centrality"
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- B[j,i]
    k = k + 1
  }
}
for (i in efforts) {
  for (j in 1:length(B[,1])) {
    r[k,1] <- "Clustering Coefficient"
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- C[j,i]
    k = k + 1
  }
}

data.stacked <- r
names(data.stacked) <- c("qtype","effort","infected")

#data.stacked = stack(data)


tx <- with(data.stacked, interaction(effort, qtype))
#analysis <- kruskal.test(data.stacked$infected, ")

# 
 data.aov = aov(data.stacked$infected~data.stacked$effort*data.stacked$qtype)
 data.aov.tk = aov(data.stacked$infected~tx)
 summary(data.aov)
 TukeyHSD(data.aov)
# 
 library(multcomp)
data <- HSD.test(data.aov.tk, "tx" ,group=T)

M <- tapply(data.stacked$infected,list(data.stacked$qtype, data.stacked$effort),mean)
sd <- tapply(data.stacked$infected,list(data.stacked$qtype, data.stacked$effort),sd)
n <- tapply(data.stacked$infected,list(data.stacked$qtype, data.stacked$effort),length)
#M <- M[,order(M[2,],decreasing=T)]

CI95 = (qnorm(0.975)  * sd / sqrt(n) )   / 3286

colors <- c("gray20",
            "gray50",
            "gray80")

bp = barplot(M/3286 ,beside=T, col=colors,
             xlab="Closure Effort (% of routes canceled)", 
             ylab="Proportion of Airports Infected",
             main="Comparison of Cancelation Strategies",
             ylim=c(0,0.5)
)
legend("topright", legend=rownames(M), fill=colors, cex=1)
plotCI(bp, as.vector(M)/3286, CI95, add=T, pch=NA)
text(x=as.vector(bp),y=as.vector(M/3286)+CI95,c("B","B","B",
                                                "D","C","A",
                                                "G","F","E",
                                                "H","J","I",
                                                "L","L","M"
                                                ), pos=3,cex=1)

abline(h=0)


# # Comparison of directed and undirected -------------------------------------
# BD <- as.data.frame(read.csv("directed/betweenness.matrix",header=F))
# BU <- as.data.frame(read.csv("undirected/betweenness.matrix",header=F))
# 
# r <- data.frame()
# k = 1
# 
# for (i in c(1,2,3,4,5)) {
#   for (j in 1:50) {
#     r[k,1] <- "betweenness"
#     r[k,2] <- as.character((i-1)*5)
#     r[k,3] <- BD[j,i]
#     r[k,4] <- "directed"
#     k = k + 1
#   }
# }
# 
# for (i in c(1,2,3,4,5)) {
#   for (j in 1:50) {
#     r[k,1] <- "betweenness"
#     r[k,2] <- as.character((i-1)*5)
#     r[k,3] <- BU[j,i]
#     r[k,4] <- "undirected"
#     k = k + 1
#   }
# }
# 
# data.stacked <- r
# names(data.stacked) <- c("qtype","effort","infected", "networktype")
# data.aov = aov(data.stacked$infected~data.stacked$networktype*data.stacked$effort)
# summary(data.aov)
# 
# # Network type is not important (P = 0.601). Effort is important (P < 0.001). 
# # There are no interations
# 
# TukeyHSD(data.aov)
# 
# tx <- with(data.stacked, interaction(effort, networktype))
# data.aov.tk = aov(data.stacked$infected~tx)
# 
# data <- HSD.test(data.aov.tk, "tx" ,group=T)
# 
# M <- tapply(data.stacked$infected,list(data.stacked$networktype, data.stacked$effort),median)
# M <- M[,order(M[2,],decreasing=T)]
# 
# CI95 = qt(1-0.05/2,9)*(sd(data.stacked$infected)/sqrt(length(data.stacked$infected)))
# 
# bp = barplot(M,beside=T, 
#              xlab="Quarantine Effort (% of airports closed)", ylab="Number of infected airports",
#              main="Comparison of betweenness-based quarantine\non different network types",
#              ylim=c(0,4500))
# legend(x=12,y=4500, legend=c("Directed","Undirected"), fill=c("gray30","gray90"), cex=0.75)
# plotCI(bp, M, CI95, add=T, pch=NA)
# text(x=as.vector(bp),y=as.vector(M)+CI95,c("A","A",
#                                            "B","B",
#                                            "C","C",
#                                            "D","D",
#                                            "D","D"), pos=3)
# abline(h=0)
# 
