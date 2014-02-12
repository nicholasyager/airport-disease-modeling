# LIST FIX

data.list <- lapply(dir(),  read.csv, header=TRUE)
main.dat <- do.call("cbind",data.list)

delete <- seq(3,134, by = 2)
main.dat <- main.dat[,-delete]



write.table(main.dat, file="look.up.table.csv",sep=",",row.names=F)
