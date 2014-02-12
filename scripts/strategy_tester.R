data = read.csv("betweenness.csv")

plot(x=data$effort, y=data$total_infected/max(data$total_infected),
     main="Effect of a betweenness-based quarantine effort\non the number of infected airports",
     xlab="Quarantine effort as proportion of all airports",
     ylab="Number of infected airports",
     pch=19)


data = read.csv("degree.csv")

plot(x=data$effort, y=data$total_infected/max(data$total_infected),
     main="Effect of a degree-based quarantine effort\non the number of infected airports",
     xlab="Quarantine effort as proportion of all airports",
     ylab="Number of infected airports",
     pch=19)
