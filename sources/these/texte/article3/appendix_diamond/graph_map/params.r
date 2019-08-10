xij <- matrix(c(
                12.26,2.18,1.64,3.19,1.48,1.09,
                2.85,23.77,2.43,3.92,2.3,1.09,
                1.33,4.29,25.2,6.71,4.83,1.56,
                6.56,6.25,11.51,25.14,15.61,4.83,
                2.18,3.13,6.95,11.45,13.89,6.38,
                0.97,1.71,2.84,5.29,6.52,9.19),
              6,6, byrow=TRUE)

cij <- matrix(c(
                169.14, 31.47, 17.76, 34.50, 15.83, 11.47,
                31.47, 274.51, 32.31, 34.86, 20.61, 11.50,
                17.76, 32.31, 224.25, 50.75, 37.52, 14.96,
                34.50, 34.86, 50.75, 75.66, 49.45, 25.08,
                15.83, 20.61, 37.52, 49.45, 61.26, 32.99,
                11.47, 11.50, 14.96, 25.08, 32.99, 54.23),
              6,6)

wi <- c(184+876,1265,1642,4857,3312,2477)*1000
wtot <- sum(wi)
Wi <- c(sum(wi[1:3]), sum(wi[4:5]), sum(wi[6]))

mij <- matrix(rep(0,6*6),6,6)
for(i in 1:6) mij[,i] <- cij[,i]*wi/wtot

Mij <- matrix(rep(0,3*3),3,3)

Mij[1,1] <-(wi[1]*(mij[1,1]+2*mij[2,1]+2*mij[3,1])+ wi[2]*(mij[2,2]+2*mij[3,2])+ wi[3]*mij[3,3])/(wi[1]+wi[2]+wi[3])
Mij[2,2] <- (wi[4]*(mij[4,4] +2*mij[5,4])+ wi[5]*mij[5,5])/(wi[4]+wi[5])
Mij[3,3] <- mij[6,6]

Mij[1,2] <- (wi[4]*(mij[1,4]+mij[2,4]+mij[3,4]) + wi[5]*(mij[1,5]+mij[2,5]+mij[3,5])) / (wi[4]+wi[5])
Mij[1,3] <- mij[1,6]+mij[2,6]+mij[3,6]
Mij[2,1] <- (wi[1]*(mij[4,1]+mij[5,1]) +wi[2]*(mij[4,2]+mij[5,2]) +wi[3]*(mij[4,3]+mij[5,3])) /(wi[1]+wi[2]+wi[3])
Mij[2,3] <-  mij[4,6]+mij[5,6]

Mij[3,1] <- (wi[1]*mij[6,1]+ wi[2]*mij[6,2] +wi[3]*mij[6,3]) /(wi[1]+wi[2]+wi[3])
Mij[3,2] <- (wi[4]*mij[6,4] +wi[5]*mij[6,5]) /  (wi[4]+wi[5])


Cij <- matrix(rep(0,3*3),3,3)
for(i in 1:3) Cij[,i] <- Mij[,i]/Wi*wtot


##############
city <- read.table("city.dat", header=TRUE)
C <- length(city[,1])

##########(min)
e <-rep(0,C)
eNS <- 0.74 ##0.8
e[city$Zone==1] <- eNS
e[city$Zone==0] <- 1
e[city$Zone==-1] <- eNS

write.table(e,file="e.dat", row.names=FALSE, col.names=FALSE)

##########
z <-rep(0,C)
z[city$Zone==1] <- 1
z[city$Zone==0] <- 1 ##(1+eNS)/2 ##on veut R0 moyen au tropiques et pas R0max
z[city$Zone==-1] <- 1

write.table(z,file="z.dat", row.names=FALSE, col.names=FALSE)

##########(tmax)
d <- rep(0,C)
d[city$Zone==1] <- 356.0  ##0
d[city$Zone==0] <- 0  ##pi/2
d[city$Zone==-1] <- 173.0  ##pi

write.table(d,file="d.dat", row.names=FALSE, col.names=FALSE)

##########
zone <-rep(0,C)
zone[city$Zone==1] <- 0
zone[city$Zone==0] <- 1
zone[city$Zone==-1] <- 2

write.table(zone,file="zone.dat", row.names=FALSE, col.names=FALSE)

##########
pai <- read.table("distrib_age.dat",header=TRUE)
pai2 <- pai[,2:4]

write.table(pai2,file="pai.dat", row.names=FALSE, col.names=FALSE)

##########
N <- as.vector(read.table("N.dat")[[1]])
Na <- matrix(rep(N,3),C,3)
Na <- pai2*N

##on fit b0 (infectiosité, q dans Wallinga) 
##Kij <- (Cij/wtot)/7
Kij <- (Cij)/7

##on calibre sur paris
KKij <- matrix(rep(0,3*3),3,3)
for(i in 1:3) KKij[,i] <- Kij[,i]*as.numeric(Na[2,])/N[2]
vp1<- eigen(KKij,symmetric=TRUE)$values

##on calibre sur mexico city
KKij <- matrix(rep(0,3*3),3,3)
for(i in 1:3) KKij[,i] <- Kij[,i]*as.numeric(Na[29,])/N[29]
vp2 <- eigen(KKij,symmetric=TRUE)$values

##tous les taux sont a exprimer en jour-1
R0 <- c(1.7, 1.7, 1.7) ##metre le premier à 0 si 2 sous type ET AUSSI I[,1]=0 Cf +bas
v <- 1/2.4 ##taux de guerison 2.4
b0 <- rep(NA,3)
b0[1] <- R0[1]*v/max(vp1)
b0[2] <- R0[2]*v/max(vp1)
b0[3] <- R0[3]*v/max(vp2)
l <- 1/1  #latence
q <- 1/(6*30.5) #taux de sortie classe Q
g <- c(1/(1*365), 1/(1*365), 1/(1*365)) ##taux de derive antigenique
tai <- c(1,1,1)
su <- c(1,1,1)

##on fabrique les fichier .dat lu par le programme C
write.table(rbind(l,v,q), file="lvq.dat", col.names=FALSE,quote=FALSE)
write.table(rbind(b0,g,tai,su), file="params.dat", col.names=FALSE,quote=FALSE)
write.table(Kij,file="Mij.dat", row.names=FALSE, col.names=FALSE)


##########
##on crée les fichiers des conditions initiales des variables d'état
##pour la co-circulation des 2 grippes saisonnieres...
##########
X <- matrix(rep(0,C*8),C,8)
E <- matrix(rep(0,C*24),C,24)
I <- matrix(rep(0,C*24),C,24)
Q <- matrix(rep(0,C*12),C,12)

################################
##modifié I1
###############################
I[,1]=20 ##I1 A MODIFIER si 2 ou 3 sous types
I[,2]=20 ##I2

X[,1] <- as.integer(Na[,1]-(I[,1]+I[,2])) ##on enleve I1+I2 de X_0
J <- cbind(X,E,I,Q)

X[,1] <- as.integer(Na[,2]-(I[,1]+I[,2]))
A <- cbind(X,E,I,Q)

X[,1] <- as.integer(Na[,3]-(I[,1]+I[,2]))
V <- cbind(X,E,I,Q)

write.table(J,file="J.dat", row.names=FALSE, col.names=FALSE)
write.table(A,file="A.dat", row.names=FALSE, col.names=FALSE)
write.table(V,file="V.dat", row.names=FALSE, col.names=FALSE)


##on calcul les R0 pour toute les villes
##R0=b/v*vp

vp3 <- rep(NA,C)
for(c in 1:C) {
  KKij <- matrix(rep(0,3*3),3,3)
  for(i in 1:3) KKij[,i] <- Kij[,i]*as.numeric(Na[c,])/N[c]
  vp3[c] <- max(eigen(KKij,symmetric=TRUE)$values)
}

r01 <- (b0[1]/v)*vp3
r02 <- (b0[2]/v)*vp3
r03 <- (b0[3]/v)*vp3

sebdb <- read.table("sebdb.dat", header=TRUE)
library("maps")

##couleur du fond des cercles
colbg <- rep('lightblue',C)
colbg[r01<=1] <- 'red'

pdf(file="r0_H3N2.pdf", width=9, height=7)
map("world")
points(sebdb[,5],sebdb[,4], pch=21, col='black', bg=colbg,cex=5*r01/max(r01))
text(sebdb[,5],sebdb[,4],round(r01,2),cex=0.7)
title("H3N2")
dev.off()

pdf(file="r0_H1N1.pdf", width=9, height=7)
map("world")
points(sebdb[,5],sebdb[,4], pch=21, col='black', bg=colbg,cex=5*r02/max(r02))
text(sebdb[,5],sebdb[,4],round(r02,2),cex=0.7)
title("H1N1")
dev.off()

pdf(file="r0_mex.pdf", width=9, height=7)
map("world")
points(sebdb[,5],sebdb[,4], pch=21, col='black', bg=colbg,cex=5*r03/max(r03))
text(sebdb[,5],sebdb[,4],round(r03,2),cex=0.7)
title("mexican flu")
dev.off()
