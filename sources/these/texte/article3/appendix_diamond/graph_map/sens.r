source("make_params.r")

##forcage saisonier (R0min/R0max dans cooper)
ee <- c(0.4, 0.74)

##############################
##R0
##on calibre sur R0 max et
##on prend R0 max dans les tropiques.
##############################

##R0 max (colonnes =sous types)
##DOIT etre une matrice sinon la fonction makeplan ne marche pas
##R0R0 <- matrix(c(
##                1.3,1.3,1.35,
##                1.5,2.8,3.5
##                )
##              ,2,3, byrow=TRUE)
##

##meme R0 pour les 3 sous types
R0R0 <- matrix(rep(seq(1.01,5,length=20),each=3), 20,3,byrow=TRUE)
##si cas 1 ou 2 sous types:
##penser a modifier make_params.r pour metre I10 ou I20 à 0
##R0R0[,1] <- 0
##R0R0[,2] <- 0


##############################
##g
##############################

##dérive antigenique (durée en année)  (colonnes =sous types)
##DOIT etre une matrice sinon la fonction makeplan ne marche pas
##gg <- matrix(c(
##                6,8,3,
##                6,8,6
##                )
##              ,2,3, byrow=TRUE)

##meme gpour les 3 sous types
gg <- matrix(rep(seq(1,10,length=20),each=3), 20,3,byrow=TRUE)


##v (TI) a rentrer en durée en jours
vv <- c(2)
##l (TE) a rentrer en durée en jours
ll <- c(1)
##durée de la classe q à rentrer en mois
qq <- c(6)

#susceptibilité des vieux
susu <- c(1)

##on genere une matrice avec en ligne une combinaison des parametres a
##faire varier et en colonne le nomre de parametres
a <- makeplan(R0R0,vv,ll,qq,gg,ee,susu)
write.table(a,file="sens.dat", row.names=FALSE, col.names=FALSE)
##on crée les fichier de parametres (fichier.numero)
for(i in 1:length(a[,1])) {
  makeparams(as.numeric(a[i,]),i)
}

