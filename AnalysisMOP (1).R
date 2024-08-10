library(fields)
library(parallel)
library(raster)
library(grDevices)
library(maptools)
library(TeachingDemos)
library(dismo)
library(biomod2)
library(sp)
library(raster)
library(rgeos)
library(rgdal)
library(usdm)
library(ENMeval)
library(foreign)
library(spocc)

rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo
#1. Crear la carpeta de la especie donde estar?n las capas de los resultados
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/")#COLOCAR AQUI LA RUTA DEL DIRECTORIO DONDE se crear? la carpeta
dir.create("Sphenarium_purpurascens")###comando para crear una carpeta con el nombre indicado entre comillas (la especie)

#########################################################################
#############Correr el analisis para cada escenario######################
#########################################################################
MOP_FOLDER <- "D:/proyectos_sigs/2._revisar_scripts/"####ruta de la carpeta donde estan todos los archivos
setwd(dir=MOP_FOLDER)
# Load functions
source("mop_v2.R")###llama a la funci?n MOP


###MOP Calibracion
#Indicar donde estan las capas de la M (presente).
m_layers_path <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/M_variables/Set_1/",
                            pattern = "*.asc$",
                            full.names = T)
#Convertir las capas de la M en un stack
mstack <- stack(m_layers_path)


###################
####Mid_Holocene###
###################
#Parte1
#Indicar donde estan las capas de la projeccion 1.
g_layers_path1 <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/G_variables/Set_1/MH_CCSM/",
                            pattern = "*.asc$",full.names = T)

#Convertir las capas de la projecci?n 1 en un stack
gstack1 <- stack(g_layers_path1)


#Correr el MOP 
mop_normalized1 <- mop(mstack,
                       gstack1,
                       percentil_prop = 0.5,
                       normalized = FALSE)  

#transformar el MOP a una capa normalizada
valor1 <- data.frame(summary(mop_normalized1))
map_fy1_mop1 <- 1-(mop_normalized1/(max(valor1[5,])))
plot(map_fy1_mop1)

##guardar el archivo
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(map_fy1_mop1, filename="Sphenarium_purpurascens_holocene1.asc", overwrite=T, suffix='names')


#Parte2
#Indicar donde est?n las capas de la projecci?n 2.
g_layers_path2 <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/G_variables/Set_1/MH_MIROC/",
                             pattern = "*.asc$",full.names = T)

#Convertir las capas de la projecci?n 1 en un stack
gstack2 <- stack(g_layers_path2)


#Correr el MOP
mop_normalized2 <- mop(mstack,
                       gstack2,
                       percentil_prop = 0.5,
                       normalized = FALSE)  

#transformar el MOP a una capa normalizada
valor2 <- data.frame(summary(mop_normalized2))
map_fy2_mop2 <- 1-(mop_normalized2/(max(valor2[5,])))

##guardar el archivo
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(map_fy2_mop2, filename="Sphenarium_purpurascens_holocene2.asc", overwrite=T, suffix='names')

#Parte3
#Indicar donde est?n las capas de la projecci?n 3.
g_layers_path3 <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/G_variables/Set_1/MH_MPI_ESM",
                             pattern = "*.asc$",full.names = T)

#Convertir las capas de la projecci?n 1 en un stack
gstack3 <- stack(g_layers_path3)

#Correr el MOP en formato normalizado
mop_normalized3 <- mop(mstack,
                       gstack3,
                       percentil_prop = 0.5,
                       normalized = FALSE)

#transformar el MOP a una capa normalizada
valor3 <- data.frame(summary(mop_normalized3))
map_fy3_mop3 <- 1-(mop_normalized3/(max(valor3[5,])))

##guardar el archivo
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(map_fy3_mop3, filename="Sphenarium_purpurascens_holocene3.asc", overwrite=T, suffix='names')

holocene_mean <- (map_fy1_mop1 + map_fy2_mop2 + map_fy3_mop3)/3
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(holocene_mean, filename="Sphenarium_purpurascens_holocene.asc", overwrite=T, suffix='names')

#1. Transformar las capas del MOP de cada a?o en mapas binarios (>0.1 NO es incertidumbre, <0.1 son ?reas de incertidumbre)
mop_holoceno_bin <- holocene_mean <= 0.1 ##selecciona como "?reas de extrapolation" aquellos sitios con valor menor a 0.1
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(mop_holoceno_bin, filename="Sphenarium_purpurascens_holocene_bin.asc", overwrite=T, suffix='names')



###############
#######LGM#####
###############
#Parte1
#Indicar donde estan las capas de la projeccion 1.
g_layers_path4 <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/G_variables/Set_1/LGM_CCSM/",
                             pattern = "*.asc$",full.names = T)

#Convertir las capas de la projecci?n 1 en un stack
gstack4 <- stack(g_layers_path4)


#Correr el MOP 
mop_normalized4 <- mop(mstack,
                       gstack4,
                       percentil_prop = 0.5,
                       normalized = FALSE)  

#transformar el MOP a una capa normalizada
valor4 <- data.frame(summary(mop_normalized4))
map_fy4_mop4 <- 1-(mop_normalized4/(max(valor4[5,])))
plot(map_fy4_mop4)

##guardar el archivo
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(map_fy4_mop4, filename="Sphenarium_purpurascens_LGM1.asc", overwrite=T, suffix='names')


#Parte2
#Indicar donde est?n las capas de la projecci?n 2.
g_layers_path5 <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/G_variables/Set_1/LGM_MIROC/",
                             pattern = "*.asc$",full.names = T)

#Convertir las capas de la projecci?n 1 en un stack
gstack5 <- stack(g_layers_path5)


#Correr el MOP
mop_normalized5 <- mop(mstack,
                       gstack5,
                       percentil_prop = 0.5,
                       normalized = FALSE)  

#transformar el MOP a una capa normalizada
valor5 <- data.frame(summary(mop_normalized5))
map_fy5_mop5 <- 1-(mop_normalized5/(max(valor5[5,])))

##guardar el archivo
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(map_fy5_mop5, filename="Sphenarium_purpurascens_LGM2.asc", overwrite=T, suffix='names')


#Parte3
#Indicar donde est?n las capas de la projecci?n 3.
g_layers_path6 <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/G_variables/Set_1/LGM_MPI_ESM",
                             pattern = "*.asc$",full.names = T)

#Convertir las capas de la projecci?n 1 en un stack
gstack6 <- stack(g_layers_path6)

#Correr el MOP en formato normalizado
mop_normalized6 <- mop(mstack,
                       gstack6,
                       percentil_prop = 0.5,
                       normalized = FALSE)

#transformar el MOP a una capa normalizada
valor6 <- data.frame(summary(mop_normalized6))
map_fy6_mop6 <- 1-(mop_normalized6/(max(valor6[5,])))

##guardar el archivo
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(map_fy6_mop6, filename="Sphenarium_purpurascens_LGM3.asc", overwrite=T, suffix='names')


LGM_mean <- (map_fy4_mop4 + map_fy5_mop5 + map_fy6_mop6)/3
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(LGM_mean, filename="Sphenarium_purpurascens_LGM.asc", overwrite=T, suffix='names')

#1. Transformar las capas del MOP de cada a?o en mapas binarios (>0.1 NO es incertidumbre, <0.1 son ?reas de incertidumbre)
mop_LGM_bin <- LGM_mean <= 0.1 ##selecciona como "?reas de extrapolation" aquellos sitios con valor menor a 0.1
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(mop_LGM_bin, filename="Sphenarium_purpurascens_LGM_bin.asc", overwrite=T, suffix='names')



###############
#######LIG#####
###############
#Parte1
#Indicar donde estan las capas de la projeccion 1.
g_layers_path7 <- list.files(path = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_purpurascens/G_variables/Set_1/LIG_CCSM/",
                             pattern = "*.asc$",full.names = T)

#Convertir las capas de la projecci?n 1 en un stack
gstack7 <- stack(g_layers_path7)


#Correr el MOP 
mop_normalized7 <- mop(mstack,
                       gstack7,
                       percentil_prop = 0.5,
                       normalized = FALSE)  

#transformar el MOP a una capa normalizada
valor7 <- data.frame(summary(mop_normalized7))
map_fy7_mop7 <- 1-(mop_normalized7/(max(valor7[5,])))
plot(map_fy7_mop7)

##guardar el archivo
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(map_fy7_mop7, filename="Sphenarium_purpurascens_LIG.asc", overwrite=T, suffix='names')

#1. Transformar las capas del MOP de cada a?o en mapas binarios (>0.1 NO es incertidumbre, <0.1 son ?reas de incertidumbre)
mop_LIG_bin <- map_fy7_mop7 <= 0.1 ##selecciona como "?reas de extrapolation" aquellos sitios con valor menor a 0.1
setwd("D:/proyectos_sigs/grasshoppers_project/MOP_final_maps/Sphenarium_purpurascens/")#DONDE Guardar la capa
writeRaster(mop_LIG_bin, filename="Sphenarium_purpurascens_LIG_bin.asc", overwrite=T, suffix='names')

