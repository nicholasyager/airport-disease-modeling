data = read.csv("betweenness.csv")

plot(data,
     main="Effect of a betweenness-based quarantine effort on the number of infected airports",
     xlab="Quarantine effort as proportion of all airports",
     ylab="Number of infected airports",
     pch=19)
