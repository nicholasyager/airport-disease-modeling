# Epidemiology of the Airways
# By Nick Yager, Matthew Taylor
# Last Update: Jan 31th, 2014
# 
#
# Purpose: To build a working network model to simulate a disease/pathogen using the 
# air travel system to propogate. 

# Libraries
library(igraph)


# Read in data files.
main.data = read.csv("airports.csv")
route.data = read.csv("routes.csv")

# Fake Vertice Fix
vertices <- 1:main.data$ID[length(main.data$ID)]
vertices <- matrix(vertices)

# Make Network Code
air.net <- graph.empty()
air.net <- add.vertices(air.net, nrow(vertices),
                        ID = as.character(vertices[,1])
                        )

edges <- matrix(c(route.data[,2], route.data[,4]), nc = 2)
edges <- edges[complete.cases(edges),] # Gets rid of NAs

air.net <- add.edges(air.net, t(edges),
                   from = edges[,1],
                   to = edges[,2]
                   )
