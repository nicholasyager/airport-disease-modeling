filelist = dir()

toremove = c()

for (i in 1:length(filelist)) {
  filesize  = file.info(filelist[i])$size
  print(filesize)
  if (filesize == 0) {
    toremove = append(toremove, i)
  }
}
filelist = filelist[-toremove]

number_of_strategies = length(filelist)

M = matrix(ncol = 21, nrow = number_of_strategies, byrow=T)

# Fill the matrix with data

for( strategy in 1:number_of_strategies ) {
  
  data = read.csv(filelist[strategy])

  for (row in 1:length(data$effort)) {
    M[strategy,row] <- data$total_infected[row]
  }
  
  
}

# Each Row is a data set for an effort. use accordingly.

# 10% effort is row 11

medians = apply(M, 2, FUN=function(v) median(v))

effort = seq(0, 100, 5)
plot(effort, medians/3307,
     main="Effect of quarantine strategies on median number of infections",
     xlab="Random quarantine effort (%)",
     ylab="Proportion of airports infected", pch=19, type="b",
     ylim=c(0,1))

abline(h=0, lty="dotted")

image(M, col=rainbow(3000))
# 
# 
#  barplot(medians/3307,
#          names=effort,
#          main="Effect of quarantine strategies on median number of infections.",
#          xlab="Degree-based quarantine effort",
#          ylab="Proportion of airports infected",
#          ylim=c(0,1))
#  abline(h=0)
