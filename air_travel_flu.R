require("igraph")
require("maps")

################### Read and prepare the route and airline data ###########################

routes = read.csv("routes.dat", header=FALSE)
airports = read.csv("airports.dat", header=FALSE)

names(airports) <- c("airport.id","name","city","country","faa","icao","latitude","longitude",
                     "alititude","timezone","dst")

names(routes) <- c("airline","airline.id","source.airport","source.id",
                   "destination.airport","destination.id","codeshare","stops","equipment")

# Limit the airports to the united states.
airports = subset(airports, country == "United States", select =  c("airport.id","name","latitude","longitude"))
airports = subset(airports, longitude < -1, select =  c("airport.id","name","latitude","longitude"))

# Limit flights to only within the United States.
routes = subset(routes, source.id %in% airports$airport.id,
                select = c("airline",
                           "airline.id",
                           "source.id",
                           "destination.airport",
                           "destination.id")
                )

routes = subset(routes, destination.id %in% airports$airport.id,
                select = c("airline",
                           "airline.id",
                           "source.id",
                           "destination.airport",
                           "destination.id")
)

################## Create the network ########################

num_airports = length(airports$airport.id)
network <- graph.data.frame( cbind(routes$source.id, routes$destination.id)  )

# Remove all verticies that do not go anywhere.
network <- delete.vertices(network, which(degree(network) < 3))
network <- delete.vertices(network, which(degree(network) < 2))
network <- delete.vertices(network, which(degree(network) < 1))

V(network)$color <- NA
E(network)$color <- NA

############### Set a random starting point

positions = c(sample(seq(1, length(V(network)),1 ), 1))
max_ticks = 1000

for (i in 1:max_ticks) {
  
  for (k in 1:length(positions)) {
    position = positions[1]
    posisitions = positions[!1]
    V(network)[position]$color <- rainbow(max_ticks)[i]
    neigh <- as.vector(V(network)[nei(position)])
    previous_position <- position
    position <- sample(neigh,1)
    E(network)[previous_position %--% position]$color <- rainbow(max_ticks)[i]
    print(position)
  }
  
}


############### Plot the network.

plot(network,
     vertex.size=5,
     vertex.label=NA,
     edge.arrow.size=0.1,
     edge.color=E(network)$color,
     vertex.color=V(network)$color
     )
