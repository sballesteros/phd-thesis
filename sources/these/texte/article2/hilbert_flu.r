#################################################
#UPCA influenza model
#################################################

#######################################
# data sentiweb
# consulté le mercredi 04 mars 2009
# Incidence (pour 100 000 habitants)
#######################################
data<-read.table("data_flu_IDF.txt")

##data filtre par fourrier lowpassfilter en enlevant 8/100 d'année
bernard <- read.table("dataflufiltre.dat")

##on interpol pour avoir plus de points
fit <- spline(seq(1,length(bernard[,1]), by=1),bernard[,1],n=100*length(bernard[,1]))
fit1 <- fit[[1]]
fit2 <- fit[[2]]
plot(fit1, fit2)
lines(seq(1,length(data[,1]), by=1),data[,1],col='red')

hilbert <- function (x)
{
    n = length(x)
    ff = fft(x)
    h = rep(0, n)
    if (n > 0 & 2 * floor(n/2) == n) {
        h[c(1, n/2 + 1)] = 1
        h[2:n/2] = 2
    }
    else {
        if (n > 0) {
            h[1] = 1
            h[2:(n + 1)/2] = 2
        }
    }
    ht = fft(ff * h, inverse = TRUE)/length(ff)
    return(ht)
}

hfit <- hilbert(fit2)
plot(Re(hfit),Im(hfit),type='l')

##attracteur UPCA sto
theo <- read.table("deter6.dat")

pdf(file="hilbert_filtered.pdf")
layout(matrix(1:2,2,1))
par(mar=c(4.0,4,1.0,0.1))
plot((Re(hfit)-mean(fit2)),Im(hfit),type='l',xlab="Re",ylab="Im")
plot(theo[1:20000,3],theo[1:20000,2],type='l',xlab="I",ylab="S")
dev.off()

