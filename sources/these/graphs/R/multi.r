postscript(file="multi.eps",width=10, height=5,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           encoding = "TeXtext.enc")

layout(matrix(1:2,1,2))
par(mar=c(4,4,2,2))
K <- seq(1,10,length=100)
plot(K,2^K, type='l', ylab=expression(2^K))

K <- seq(1,40,length=100)
plot(K,2^K, type='l', ylab=expression(2^K))

dev.off()
