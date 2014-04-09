# Map.R
# Nicholas A. Yager and Matthew Taylor.
#
# An early attempt to plot maps of the network overlayed onto a projection of the world

library("maps")

airports = read.csv("data/airports.dat")
data = read.csv("betweenness.csv")
colorValues = data$betweenness*10000

ramp <- colorRampPalette(c("blue","red"))
colors = ramp(max(colorValues))

data = data[order(data$betweenness),]

map("world", xlim=c(-180,180), ylim=c(-90,90), mar=c(0,0,0,0))

points(x=data$lon, y=data$lat, pch=16, cex=0.4, col=colors[colorValues])
