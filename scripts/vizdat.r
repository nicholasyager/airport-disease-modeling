# Purpose: Visualize the data for random vaccination efforts, and be able to output the variance,
# median, mean, and eventually build an ANOVA. All data to be sorted by vaccination effort.
# I: csv Data sets
# O: Graphics, Look-up Matrix, and Stats.

library(scales)

# Attach Data
main.data = read.csv("betweenness.table.csv")

# Add in points.
plot(main.data[,1],main.data[,2],
     col = alpha("black",0.4), pch = 16,
     ylab = "Number of Airport Closures",
     xlab = "Quarantine Effort (% Closures)"
     )
for (i in 1:(ncol(main.data) - 1)){
  points(main.data[,1], main.data[,i + 1],
         col = alpha("black",0.4), pch = 16
         )
}

median.matrix <- main.data[,-1]
median.matrix <- t(median.matrix)
med.vals <- F
mean.vals <- F
for (i in 1:(ncol(median.matrix))){
  med.vals[i-1] <- median(median.matrix[,i])
  mean.vals[i-1] <- mean(median.matrix[,i])
}

med.vals[101] <- 0
mean.vals[101] <- 0

lines(main.data[,1],mean.vals, col = "red", lwd = 1.5)
lines(main.data[,1],med.vals, col = "blue", lwd = 1.5)

legend("right", c("Mean", "Median"), col = c("red", "blue"),
       lty = 1, lwd = 1.5
       )

boxplot(t(main.data))

lines(1:nrow(main.data),mean.vals, col = "red", lwd = 1.5)
lines(1:nrow(main.data),med.vals, col = "blue", lwd = 1.5)

legend("right", c("Mean", "Median"), col = c("red", "blue"),
       lty = 1, lwd = 1.5
)

