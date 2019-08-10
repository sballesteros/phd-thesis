sigma <- 0.1
traj.gupta <- read.table(paste("gupta/traj.",sigma,sep=""))
traj.sbri <- read.table(paste("sb_ri/traj.",1-sigma,sep=""))
traj.sbrs <- read.table(paste("sb_rs/traj.",1-sigma,sep=""))
traj.hbri <- read.table(paste("hb_ri/traj.",sigma,sep=""))
traj.hbrs <- read.table(paste("hb_rs/traj.",sigma,sep=""))

maxI2 <- max(c(max(traj.gupta[,3]),max(traj.sbrs[,3]),max(traj.sbri[,3]),max(traj.hbrs[,3]), max(traj.hbri[,3])))
minI2 <- min(c(min(traj.gupta[,3]),min(traj.sbrs[,3]),min(traj.sbri[,3]),min(traj.hbrs[,3]), min(traj.hbri[,3])))

pdf(file="compare1.pdf")
plot(traj.gupta[,1],log10(traj.gupta[,2]),type='l',lty=2, col='black',
     xlim=c(0,10),
     ylim=log10(c(minI2,maxI2)),
     xlab="time (years)",
     ylab="log10(I)")
     mtext("SBRS, HBRS, HBRI, Gupta, ",  side = 3, line = 1, cex = 2, at=4, col = 'black')
lines(traj.gupta[,1],log10(traj.gupta[,3]), col='black')
lines(traj.sbrs[,1],log10(traj.sbrs[,2]),lty=2, col='black')
lines(traj.sbrs[,1],log10(traj.sbrs[,3]), col='black')
lines(traj.hbri[,1],log10(traj.hbri[,2]),lty=2, col='black')
lines(traj.hbri[,1],log10(traj.hbri[,3]), col='black')
lines(traj.hbrs[,1],log10(traj.hbrs[,2]),lty=2, col='black')
lines(traj.hbrs[,1],log10(traj.hbrs[,3]), col='black')
legend("topright",c("redident" ,"mutant (10 % im. escape)"), lty=c(2,1))
dev.off()

pdf(file="compare2.pdf")
plot(traj.gupta[,1],log10(traj.gupta[,2]),type='l',lty=2, col='black',
     xlim=c(0,10),
     ylim=log10(c(minI2,maxI2)),
     xlab="time (years)",
     ylab="log10(I)")
     mtext("SBRS, HBRS, HBRI, Gupta, ",  side = 3, line = 1, cex = 2, at=4, col = 'black')
     mtext("SBRI",  side = 3, line = 1, cex = 2, at=9, col = 'red')
lines(traj.gupta[,1],log10(traj.gupta[,3]), col='black')
lines(traj.sbrs[,1],log10(traj.sbrs[,2]),lty=2, col='black')
lines(traj.sbrs[,1],log10(traj.sbrs[,3]), col='black')
lines(traj.hbri[,1],log10(traj.hbri[,2]),lty=2, col='black')
lines(traj.hbri[,1],log10(traj.hbri[,3]), col='black')
lines(traj.hbrs[,1],log10(traj.hbrs[,2]),lty=2, col='black')
lines(traj.hbrs[,1],log10(traj.hbrs[,3]), col='black')
lines(traj.sbri[,1],log10(traj.sbri[,2]),lty=2, col='red',lwd=2)
lines(traj.sbri[,1],log10(traj.sbri[,3]), col='red',lwd=2)
legend("topright",c("redident" ,"mutant (10 % im. escape)"), lty=c(2,1))
dev.off()
