###cartographie...
library("maps")
data(world.cities)
trans <- read.table("/home/seb/Bureau/pandemics/drift/data/transport2000.dat", header=TRUE)

##on remplace Hong kong par Shenzhen car pas de data dans Rmap
city <- read.table("/home/seb/Bureau/pandemics/drift/data/city_renamed.dat",header=TRUE)
##city <- city[-10,]
country <- read.table("/home/seb/Bureau/pandemics/drift/data/coutry.dat")[[1]]

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

flux <- data.frame(from=trans[trans$flow !=0,1],to=trans[trans$flow !=0,2], flux=trans[trans$flow !=0,3]/max(trans[trans$flow !=0,3]))

#pdf(file="map.pdf", width=20,height=15)

postscript(file="graph/iata.eps",width=8.267, height=7,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           encoding = "TeXtext.enc")

layout(matrix(1:2,2,1))


postscript(file="graph_annexe/iata.eps",width=8.267, height=4,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           encoding = "TeXtext.enc")

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


dev.off()



##depart depuis mexico en vert
for(i in 442:463) {
  x <- c(sebdb[sebdb$ind==trans2$from2[i],5], sebdb[sebdb$ind==trans2$to2[i],5])
  y <- c(sebdb[sebdb$ind==trans2$from2[i],4], sebdb[sebdb$ind==trans2$to2[i],4])
  lines(x,y, lwd=10000*prob[i], col='green')
}

pdf(file="metapop.pdf")

I1 <- read.table("../metapop/src/J1.0.dat", header=TRUE)
I2 <- read.table("../metapop/src/J2.0.dat", header=TRUE)
t <- seq(1,length(I1[,1]),by=1)/52

tmax <- 25
city <- 2
par(mar=c(4,4,2,2))
plot(t[t<tmax], I1[t<tmax,city]+I2[t<tmax,city], type='l', xlab="time (years)", ylab="I1,I2,I1+I2", lwd=2 )
lines(t[t<tmax], I1[t<tmax,city], col='red' )
lines(t[t<tmax], I2[t<tmax,city], col='blue' )


dev.off()
