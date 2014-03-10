library(plotrix)
library(agricolae)

R <- as.data.frame(read.csv("random.matrix",header=F))
#D <- as.data.frame(read.csv("degree.matrix",header=F))
B <- as.data.frame(read.csv("betweenness.matrix",header=F))

r <- data.frame()
k = 1

efforts = c(1,3,5,7,9)

for (i in efforts) {
  for (j in 1:length(R[,1])) {
    r[k,1] <- 1
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- R[j,i]
    k = k + 1
  }
}
# for (i in c(1,3,5,7,8)) {
#   for (j in 1:50) {
#     r[k,1] <- 2
#     r[k,2] <- as.character((i-1)*5)
#     r[k,3] <- D[j,i]
#     k = k + 1
#   }
# }
for (i in efforts) {
  for (j in 1:length(B[,1])) {
    r[k,1] <- 3
    r[k,2] <- as.character((i-1)*5)
    r[k,3] <- B[j,i]
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

M <- tapply(data.stacked$infected,list(data.stacked$qtype, data.stacked$effort),median)
M <- M[,order(M[2,],decreasing=T)]

CI95 = qt(1-0.05/2,9)*(sd(data.stacked$infected)/sqrt(length(data.stacked$infected))) / 3287

bp = barplot(M/3287 ,beside=T,
             xlab="Quarantine Effort (% of edges closed)", 
             ylab="Proportion of infected airports",
             main="Comparison of Quarantine Strategies",
             ylim=c(0,1)
)
legend(x=12,y=1, legend=c("Random","Betweenness"), fill=c("gray30","gray90"), cex=0.75)
plotCI(bp, M/3287, CI95, add=T, pch=NA)
text(x=as.vector(bp),y=as.vector(M/3287)+CI95,c("A","A",
                                           "B","C",
                                           "D","E",
                                           "F","G",
                                           "H","I"
                                           ), pos=3,cex=0.8)

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
