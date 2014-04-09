library("igraph")

edgelist <- read.csv("edgelist.csv")
airport_data <- read.csv("airports.dat",header=F)
  
relations <- data.frame(as.character(edgelist$IATA_From), as.character(edgelist$IATA_To))
names <- append(as.character(edgelist$IATA_From), as.character(edgelist$IATA_To))
duplicates <- duplicated(names)
names <- data.frame(name=names[!duplicates])

g <- graph.data.frame(relations, directed=T, vertices=names)



which(degree(g) < 2)

# Degree distribution for all airports
hist(log(degree(g)),
     main="Distribution of Network Degree",
     xlab="ln(Degree)",
     col="gray70")

# Degree distribution for ATL
ATL_data <- log(degree(g,V(g)[neighbors(g,'ATL')]))
his <- hist(log(degree(g,V(g)[neighbors(g,'ATL')])),
     main="Degree distribution of ATL's Neighbors",
     col="firebrick4",
     xlab="ln(Degree)",
     xlim=c(min(ATL_data),max(ATL_data)))

his2 <- hist(log(degree(g,V(g)[neighbors(g,'ROC')])),
     col="cadetblue4",
     main="Degree distribution of ROC's Neighbors",
     xlab="ln(Degree)",
     xlim=c(min(ATL_data),max(ATL_data)), 
     ylim=c(0, max(his$counts)),
     breaks=his$breaks)

