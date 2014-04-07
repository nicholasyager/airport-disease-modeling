# Delay Aggregator

library("agricolae")
library("plotrix")

# Fetch list of directories
dirs <- dir()

# Remove the unwanted directories
internationals <- dir(pattern="international") # remove these
to_remove <- match(internationals,dirs)
dirs <- dirs[-to_remove]

# Load the list of desired factors
factors <- c("betweenness","cluster")

data <- matrix(ncol=4)
# Iterate through the folders
for(dir_i in 1:length(dirs)){
  
  # Recurse into the factor folders
  for(factor_i in 1:length(factors)) {
   
    file_dir = paste("./",dirs[dir_i],"/",factors[factor_i],sep="")
    files <- dir(file_dir)
    
    # Read each file and place it in the master data.frame
    for(file_i in 1:length(files)){
      file_path <- paste(file_dir,"/",files[file_i],sep="")
      file_data <- read.csv(file_path)
      
      for (row_i in 1:length(file_data[,1])) {
        new_data <- c( as.numeric(dir_i - 1), factors[factor_i], as.numeric(file_data$effort[row_i]), file_data$total_infected[row_i])
        print(new_data)
        data <- rbind(data,new_data)
      } 
    }
  }
}

data <- data[-1,]
row.names(data) <- NULL
data.frame <- as.data.frame(data)
names(data.frame) <- c("delay","qtype","effort","infected")

test <- subset(data.frame, effort == "0.31")
row.names(test) <- NULL
test$infected <- as.numeric(test$infected)
test$effort <- test$effort

tx <- with(test, interaction(qtype,delay))
data.aov <- aov(test$infected~tx)

HSD <- HSD.test(data.aov, "tx" ,group=T)

M <- tapply(test$infected,list(test$qtype, test$delay),mean)
sd <- tapply(test$infected,list(test$qtype, test$delay),sd)
n <- tapply(test$infected,list(test$qtype, test$delay),length)

CI95 = (qnorm(0.975)  * sd / sqrt(n) )   / 3286

colors <- c("gray40",
            "gray80")

bp = barplot(M/3286 ,beside=T, col=colors,
             xlab="Delay (Days)", 
             ylab="Proportion of Airports Infected",
             main="Comparison of Strategies with Delays at 30% Effort",
             ylim=c(0,0.5)
)
legend("topright", legend=c("Betweenness Centrality","Clustering Coefficient"), fill=colors, cex=0.8)
plotCI(bp, as.vector(M)/3286, CI95, add=T, pch=NA)
text(x=as.vector(bp),y=as.vector(M/3286)+CI95,c("H","FG",
                                                "GH","FG",
                                                "H","EF",
                                                "H","E",
                                                "D","E",
                                                "D","B",
                                                "A","C"
), pos=3,cex=0.8)
par(xpd=F)
abline(h=0)
