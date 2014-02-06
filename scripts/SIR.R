data = read.csv("data.csv")

max = max(data$s)

plot(data$s/max, 
     main="Simplistic infection simulation with an airport-based network model",
     xlab="Time (Tick)",
     ylab="Proportion of Airports",
     type="l",
     col="blue")
lines(data$i/max, type="l",
      col="green",
      lty="dashed")
lines(data$r/max, type="l",
      col="red",
      lty="dotted")
legend(x=11.5, y=0.6, legend=c("Succeptable","Infected","Recovered"), 
       lty=c("solid","dashed","dotted"),
       col=c("blue","green","red"))
