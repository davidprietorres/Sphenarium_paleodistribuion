library(dismo)
library(ade4)
library(sp)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(usdm)
library(foreign)
library(spocc)
library(corrplot)
library(usdm)
library(XML)
library(dplyr)
library(raster)
library(dismo) 
library(ecospat)   
library(phyloclim)


rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo

##1.Read the .CSV file to the occurrence point for the species
species_1<- read.csv("D:/proyectos_sigs/salomon_diversiﬁcation of grasshoppers/Sampling records/Sphenarium variabile.csv", header = T, sep = ",")#occurrence for the especie 1
names(species_1)
species_1$lon <- as.numeric(species_1$lon)
species_1$lat <- as.numeric(species_1$lat)

species_2 <- read.csv("D:/proyectos_sigs/salomon_diversiﬁcation of grasshoppers/Sampling records/Sphenarium rugosum.csv", header = T, sep = ",")#occurrence for the especie 2
species_2$lon <- as.numeric(species_2$lon)
species_2$lat <- as.numeric(species_2$lat)

##2.Read the environmental variables information for each species
##Climate information for all study area
setwd("D:/proyectos_sigs/salomon_diversiﬁcation of grasshoppers/capas_climaticas/presente2.1/") 
variables <- list.files(".",pattern = "*.asc$",full.names = T)###crea el stack 
varclim<- stack(variables)

#M specie1
map_m_spp1 <- readOGR("D:/proyectos_sigs/salomon_diversiﬁcation of grasshoppers/M_analyses/Simulations/Sphenarium_variabile.shp")
plot(map_m_spp1)

#M specie2
map_m_spp2 <- readOGR("D:/proyectos_sigs/salomon_diversiﬁcation of grasshoppers/M_analyses/Simulations/Sphenarium_rugosum.shp")
plot(map_m_spp2)

#Clip the environmental variables to the M extent.
##species1
CROPPED_m_spp1 <- raster::crop(varclim, extent(map_m_spp1))
varclim_sp1 <- mask(CROPPED_m_spp1, map_m_spp1)
plot(varclim_sp1)

##species2
CROPPED_m_spp2 <- raster::crop(varclim, extent(map_m_spp2))
varclim_sp2 <- mask(CROPPED_m_spp2, map_m_spp2)
plot(varclim_sp2)


##3.To create a data frame with the available climate conditions for each M area
clim_punto_sp1 <- rasterToPoints(varclim_sp1[[1]], fun=NULL, spatial=TRUE)##puntos background de la especie 1
clim_punto_sp2 <- rasterToPoints(varclim_sp2[[1]], fun=NULL, spatial=TRUE)##puntos background de la especie 2


##4. To extract the climate information for each coordinate
##species 1
clima_species_1 <- raster::extract(varclim_sp1, clim_punto_sp1)
clima_species_1 <- data.frame(coordinates(clim_punto_sp1),clima_species_1)
clima_species_11 <- subset(clima_species_1, !is.na(bio_1) & !is.na(bio_2) & !is.na(bio_3) & !is.na(bio_8) & !is.na(bio_11) & !is.na(bio_12) & !is.na(bio_14)& !is.na(bio_15) & !is.na(bio_17))
names (clima_species_1)[1] = "lon"
names (clima_species_1)[2] = "lat"

##species 2
clima_species_2 <- raster::extract(varclim_sp2, clim_punto_sp2)
clima_species_2 <- data.frame(coordinates(clim_punto_sp2),clima_species_2)
clima_species_22 <- subset(clima_species_2, !is.na(bio_1) & !is.na(bio_2) & !is.na(bio_3) & !is.na(bio_8) & !is.na(bio_11) & !is.na(bio_12) & !is.na(bio_14)& !is.na(bio_15) & !is.na(bio_17))
names (clima_species_2)[1] = "lon"
names (clima_species_2)[2] = "lat"

#5. To selecte the data presence (only the coordinates for each species)
occ.sp1 <- species_1[2:3]
occ.sp2 <- species_2[2:3]

#6.Integrating the information of species occurrence and environmental conditions for each species in only file
occ_sp1 <- na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clima_species_1,
                                         colvarxy=1:2,colvar="all",resolution= 0.008333))

occ_sp2 <- na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                         colspkept=1:2,dfvar=clima_species_2, 
                                         colvarxy=1:2,colvar="all",resolution=0.008333))


#7. To include a new row with the information for the category of data: occurrence (1) vs. background (0)
#specie1
occ_sp1 <- cbind(occ_sp1,species_occ = 1)###agrega una columna con valor "1" para los datos de presencia de la sp1.
clima_species_1 <- cbind(clima_species_1,species_occ = 0)###agrega una columna con valor "0" para los datos de background (M) de la sp1.)

names(clima_species_1)[1] = "lon" ##cambia el nombre de la columna longitud para que coincida con el nombre del archivo de presencias
names(clima_species_1)[2] = "lat" ##cambia el nombre de la columna latitud para que coincida con el nombre del archivo de presencias
data_species1 = rbind(occ_sp1, clima_species_1)##unir los datos de presencia y background para la species 1


#specie2
occ_sp2 <- cbind(occ_sp2,species_occ = 1)###agrega una columna con valor "1" para los datos de presencia de la sp1.
clima_species_2 <- cbind(clima_species_2,species_occ = 0)###agrega una columna con valor "0" para los datos de background (M) de la sp1.)

names(clima_species_2)[1] = "lon" ##cambia el nombre de la columna longitud para que coincida con el nombre del archivo de presencias
names(clima_species_2)[2] = "lat" ##cambia el nombre de la columna latitud para que coincida con el nombre del archivo de presencias
data_species2 = rbind(occ_sp2, clima_species_2)##unir los datos de presencia y background para la species 2


##8. To perfom the PCA with the information for both species
pca.env <- dudi.pca(rbind(data_species1,data_species2)[,3:11],center = T, scale = T, scannf=F,nf=2)
summary(pca.env)

windows()
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)###grafica el PCA obtenido mostrando la contribuci?n de los dos componentes principales

scores.globclim <- pca.env$li

scores.sp1 <- suprow(pca.env,data_species1[which(data_species1[,12]==1),3:11])$li#PCA scores for the specie 1
scores.sp2 <- suprow(pca.env,data_species2[which(data_species2[,12]==1),3:11])$li#PCA scores for the specie 2.

scores.clim.sp1 <- suprow(pca.env,data_species1[which(data_species1[,12]==0),3:11])$li#PCA scores for the background area of species 1
scores.clim.sp2 <- suprow(pca.env,data_species2[which(data_species2[,12]==0),3:11])$li#PCA scores for the background area of species 1


#9.To stimate the graph of density of species' occurrence
#Species1
grid.clim.sp1 <- ecospat.grid.clim.dyn(glob= scores.globclim,
                                       glob1= scores.clim.sp1,
                                       sp=scores.sp1, R=300,
                                       th.sp=0)

#Species2
grid.clim.sp2 <- ecospat.grid.clim.dyn(glob=scores.clim.sp2,
                                              glob1=scores.clim.sp2,
                                              sp=scores.sp2, R=300,
                                              th.sp=0)


###Calculate Niche Overlap with ecospat.niche.overlap()Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.sp1, grid.clim.sp2, cor= T)
D.overlap###este es el valor de sobrelape observado (valor empir?co)... es decir la flecha a comparar en las gr?ficas! Este valor es el que se reporta en el Excel.


##para ver los gr?ficos de nicho!
windows() # en mac o windows() en PC
###par(mfrow=c(1,3))

##ver las graficas de los nichos en la misma figura.
ecospat.plot.niche.dyn(grid.clim.sp1, grid.clim.sp2, quant=0.25,colZ1="Green", colZ2 = "Red", 
                       title= "Niche Overlap", name.axis1="PC1",name.axis2="PC2")


##realizar los analisis de equivalencia de nicho
equiv.test <- ecospat.niche.equivalency.test(grid.clim.sp1, grid.clim.sp2,
                                             rep=1000)



#11. Realizar los an?lisis de similitud de nichos ("Background similiarity test")
#Comparaci?n sp1 vs. sp2.
sim.test_sp1 <- ecospat.niche.similarity.test(grid.clim.sp1, grid.clim.sp2,
                                            rep=1000, 
                                            rand.type=2)##esta hecho para hacer mil replicas comparando los datos 
sim.test_sp2 <- ecospat.niche.similarity.test(grid.clim.sp2, grid.clim.sp1,
                                            rep=1000, 
                                            rand.type=2)##esta hecho para hacer mil replicas comparando los datos

windows() # en mac o windows() en PC
par(mfrow=c(1,3))
ecospat.plot.overlap.test(equiv.test,"D","Equivalence - S. variabile vs. S. rugosum")##Cambiar los nombres de los titulos. El valor de "P" en la gr?fica es lo que se reporta en el Excel.
ecospat.plot.overlap.test(sim.test_sp1,"D","Similarity - S. variabile vs. S. rugosum")##Cambiar los nombres de los titulos. El valor de "P" en la gr?fica es lo que se reporta en el Excel.
ecospat.plot.overlap.test(sim.test_sp2,"D", "Similarity - S. rugosum vs. S. variabile")##Cambiar los nombres de los titulos. El valor de "P" en la gr?fica es lo que se reporta en el Excel.

###FIN