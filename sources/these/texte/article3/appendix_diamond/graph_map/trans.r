city <- read.table("city.dat", header=TRUE)
N <- length(city[,1])
trans <- read.table("transport2000.dat", header=TRUE)
##N <- 4
##trans <- read.table("transport_test.dat", header=TRUE)
attach(trans)
##flux par jour

##on reconstruit les liens complet en supposant symetrie des fluxs
trans2 <- data.frame(from2=trans[(from==1),1], to2=trans[(from==1),2], flow2=trans[(from==1),3])

for(i in 2:N){
  transtemp.b <- data.frame(from2=trans[(to==i),2],
                            to2=trans[(to==i),1],
                            flow2=trans[(to==i),3])
  transtemp.f <- data.frame(from2=trans[(from==i),1],
                            to2=trans[(from==i),2],
                            flow2=trans[(from==i),3])
  trans2 <- rbind(trans2,transtemp.b,transtemp.f)
}

trans2 <- trans2[trans2[,3]!=0.0,]

##nombre de villes qu'on peut atteindre depuis chacques villes
nb <- table(as.factor(trans2[,1]))

write.table(as.vector(nb),file="depart.dat",
            row.names = FALSE, col.names = FALSE)

##indices de départ
indd <- seq(1, nb[1],by=1)
for(i in 2: length(nb)){
  tmp <- seq(1, nb[i],by=1)
  indd <- c(indd,tmp)
}

trans2 <- data.frame(trans2,indd=indd)

##indices d'arrivés
inda <- trans2[,4]
for(i in 1:length(inda)){
  inda[i] <- trans2[(trans2$from2==trans2[i,2]) & (trans2$to2==trans2[i,1]), 4]
}

trans2 <- data.frame(trans2,inda=inda)

reach <- data.frame(arrive=trans2[,1], depart=trans2[,2], reaction=trans2[,5])

write.table(reach,file="reach.dat",
            row.names = FALSE, col.names = FALSE)

##on transform les flux en probas de passages en divisant par la taille de la pop de départ
city2 <- data.frame(city=seq(1,length(city[,1]),by=1), pop=city$Populationsize*1000)

prob <- trans2$flow2
for(i in 1:length(prob)){
 prob[i] <- prob[i]/city2$pop[city2$city==trans2$from2[i]]
}

write.table(prob,file="prob_trans.dat",
            row.names = FALSE, col.names = FALSE)

write.table(as.integer(city2$pop),file="N.dat",
            row.names = FALSE, col.names = FALSE)


###cartographie...
library("maps")
data(world.cities)

##on remplace Hong kong par Shenzhen car pas de data dans Rmap
city <- read.table("city_renamed.dat",header=TRUE)
##city <- city[-10,]
country <- read.table("coutry.dat")[[1]]

city <- data.frame(City=city[,1], Country=country, Populationsize=city[,2]*1000, Zone=city[,3])

##A list with 6 components, namely "name", "country.etc", "pop", "lat", "long", and "capital", containing the city name, the country name, approximate population (as at January 2006), latitude,longitude and capital status indication (0 for non-capital, 1 for capital, 2 for China Municipalities,and 3 for China Provincial capitals)
lat <- rep(NA,length(city[,1]))
long <- rep(NA,length(city[,1]))
capital <- rep(NA,length(city[,1]))

for(i in 1:length(city[,1])) {
lat[i] <- world.cities[(world.cities$name==as.character(city[i,1])) & (world.cities$country.etc==as.character(city[i,2])),4]
long[i] <- world.cities[(world.cities$name==as.character(city[i,1])) & (world.cities$country.etc==as.character(city[i,2])),5]
capital[i] <- world.cities[(world.cities$name==as.character(city[i,1])) & (world.cities$country.etc==as.character(city[i,2])),6]
##print(world.cities[(world.cities$name==as.character(city[i,1])) & (world.cities$country.etc==as.character(city[i,2])),1])
}

sebdb <- data.frame(name=as.character(city[,1]), country.etc=city[,2], pop=city[,3], lat=lat, long=long, zone=city[,4], ind=seq(1,length(city[,1]),by=1))

write.table(sebdb,file="sebdb.dat", row.names = FALSE)

flux <- data.frame(from=trans[trans$flow !=0,1],to=trans[trans$flow !=0,2], flux=trans[trans$flow !=0,3]/max(trans[trans$flow !=0,3]))

postscript(file="mapHK.eps",width=10, height=7.5,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           encoding = "TeXtext.enc")

par(mar=c(1,1,1,1))
map("world")
##map.cities(x=sebdb,label=FALSE,col='red',cex=2,fill=TRUE)
points(sebdb[,5],sebdb[,4], pch=21, bg='red',cex=8*sebdb[,3]/max(sebdb[,3]))

for(i in 1:length(flux[,1])) {
  x <- c(sebdb[sebdb$ind==flux$from[i],5], sebdb[sebdb$ind==flux$to[i],5])
  y <- c(sebdb[sebdb$ind==flux$from[i],4], sebdb[sebdb$ind==flux$to[i],4])
  lines(x,y, lwd=2*flux[i,3])
}

##depart depuis HK en vert
for(i in 179:200) {
  x <- c(sebdb[sebdb$ind==trans2$from2[i],5], sebdb[sebdb$ind==trans2$to2[i],5])
  y <- c(sebdb[sebdb$ind==trans2$from2[i],4], sebdb[sebdb$ind==trans2$to2[i],4])
  lines(x,y, lwd=10000*prob[i], col='green')
}


dev.off()



##pres these


pdf(file="map2.pdf", width=10,height=5)
par(mar=c(0,0,0,0))
cz <- rep(NA, length(sebdb$zone))
cz[sebdb$zone==1] <- 'red'
cz[sebdb$zone==0] <- 'green'
cz[sebdb$zone==-1] <- 'blue'

par(mar=c(0,0,0,0))
map("world")
##map.cities(x=sebdb,label=FALSE,col='red',cex=2,fill=TRUE)
points(sebdb[,5],sebdb[,4], pch=21, bg=cz,cex=8*sebdb[,3]/max(sebdb[,3]))

for(i in 1:length(flux[,1])) {
  x <- c(sebdb[sebdb$ind==flux$from[i],5], sebdb[sebdb$ind==flux$to[i],5])
  y <- c(sebdb[sebdb$ind==flux$from[i],4], sebdb[sebdb$ind==flux$to[i],4])
  lines(x,y, lwd=2*flux[i,3])
}


##depart depuis HK en vert
for(i in 179:200) {
  x <- c(sebdb[sebdb$ind==trans2$from2[i],5], sebdb[sebdb$ind==trans2$to2[i],5])
  y <- c(sebdb[sebdb$ind==trans2$from2[i],4], sebdb[sebdb$ind==trans2$to2[i],4])
  lines(x,y, lwd=10000*prob[i], col='cyan')
}

dev.off()
