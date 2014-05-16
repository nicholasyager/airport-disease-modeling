require("igraph")

trim.trailing <- function (x) sub("\\s+$","",x)

# FAA Data
data <- read.csv("TFMSC-Report-84716.csv", stringsAsFactors=F, strip.white=T)
data$IATA <- trim(data$IATA)

# Network Data
edgelist <- read.csv("edgelist.csv")
edges <- matrix(c(as.character(edgelist$IATA_From), as.character(edgelist$IATA_To)), ncol=2)
G <- graph.edgelist(edges, directed=T)

# Find degrees
degrees <- degree(G)

degrees <- data.frame(IATA=as.character(names(degrees)), degree=degrees, stringsAsFactors=F)
row.names(degrees) <- NULL

# Cross-reference

output <- matrix(nrow=1, ncol=5)

# For each airport in the US dataset
for(i in 1:length(data$IATA)) {
  
  # Look for a correspoding degree
  for (j in 1:length(degrees$IATA)) {
    
    if (data$IATA[i] == degrees$IATA[j]) {
      if (i == 1) {
        output[i,] <- c(data$IATA[i], degrees$degree[j], data$Depart.AC.AT[i], data$Arrival.AC.AT[i], data$Total[i])
      } else {
        output <- rbind(output, c(data$IATA[i], degrees$degree[j], data$Depart.AC.AT[i], data$Arrival.AC.AT[i], data$Total[i]))
      }
      print(data$IATA[i])
      break
    }
    
  }
  
}

hist(as.numeric(output[,2]))
hist(as.numeric(output[,5]))

plot(as.numeric(output[,2]),as.numeric(output[,5]),
     xlab="Airport Degree",
     ylab="Average number of commercial operations per day",
     main="The relationship between airport degree and the average\n number of daily operations is linearly correlated",
     pch=19, cex=0.75
     )
x <- as.numeric(output[,2])
y <- as.numeric(output[,5])


mod <- lm(y ~ x )
summary <- summary(mod)

abline(mod, col="red")
text(325, 20000,paste("y =",mod$coefficients[2],"x +",mod$coefficients[1],"\rR-squared =",summary$adj.r.squared))
