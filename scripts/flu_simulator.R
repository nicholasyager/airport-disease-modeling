# Simulator.R
#
# Matthew Taylor and Nicholas A. Yager

library("igraph")

# Load the airport data.
airports = read.csv("data/airports-clean.csv")
names(airports) = c("table_id","id","name","city","country","FAA","IATA","lat","long","alt","timezone","dst")
routes = read.csv("data/routes-clean.csv")
names(routes) = c("table_id","airline","airline_id","source_IATA","source_id","dest_IATA","dest_id","codeshare","stops","equipment")

# Clear the tables of malformed results
routes$source_id <- as.numeric(as.character(routes$source_id))
routes$dest_id <- as.numeric(as.character(routes$dest_id))

routes = routes[complete.cases(routes[,c(5,7)]),]

# Add edges and vertices
edgelist <- cbind(routes$source_id, routes$dest_id)
network = graph(edgelist,directed=T)

degree = centralization.degree(network, mode="in")
for (i in 1:length(degree$res)) {
  if (degree$res[i] < 1) {
    print(i)
    v = V(network)[i]
    delete.vertices(network,v)
    
  }
}

start = 3682
S = 0
I = 1
R = 0


V(network)$status <- "S"
V(network)[start]$status <- "I"


for (i in 1:length(V(network))) {
  if (V(network)[i]$status == "I")
  {
    neighbors = neighbors(network,V(network)[i])
    print(i)
    print(neighbors)
  }
}

lo <- layout.auto(network)
plot(network, layout=lo)