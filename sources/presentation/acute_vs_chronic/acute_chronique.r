traja <- read.table("traj_acute.dat")
trajc <- read.table("traj_chronique.dat")



pdf(file="acute0.pdf", width=9,height=6)
par(mar=c(4,4,0,0))
plot(trajc[,1]-20,trajc[,3], type='l',
     xlab="time (years)", ylab="I",
     xlim=c(20,60), ylim=c(0,12000))
legend("topright", c("chronic: R0=5, duration of infection=1 year"), lty=1)
dev.off()

pdf(file="acute1.pdf", width=9,height=6)
par(mar=c(4,4,0,0))
plot(trajc[,1]-20,trajc[,3], type='l',
     xlab="time (years)", ylab="I",
     xlim=c(20,60), ylim=c(0,12000))
lines(traja[,1],traja[,3], col='red')
legend("topright", c("chronic: R0=5, duration of infection=1 year", "acute: R0=5, duration of infection=7 days"), lty=1, col=c('black','red'))
dev.off()
