library(rgbif)
library(TeachingDemos)
library(dismo)
library(biomod2)
library(sp)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(usdm)
library(ENMeval)
library(foreign)
library(spocc)
library(corrplot)
library(usdm)
library(XML)
library(dplyr)
library(reshape)
library(kuenm)

rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo
##################################################
##CUT THE BIOCLIMATE VARIABLES INTO "M" EXTENT####
##################################################
#PARTE I: Create the folders and directories
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/")#
dir.create("Sphenarium_zapotecum")###comando para crear una carpeta con el nombre indicado entre comillas

setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/")
dir.create("G_variables")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("M_variables")###comando para crear una carpeta con el nombre indicado entre comillas

#3.Crear las carpetas para los dos "Set" diferentes de variables del presente (Set1)
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/M_variables/")#
dir.create("Set_1")###comando para crear una carpeta con el nombre indicado entre comillas

#4.Crear las carpetas para los dos "Set" diferentes de variables del futuro (Set1)
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/")#
dir.create("Set_1")###comando para crear una carpeta con el nombre indicado entre comillas

#5.Crear las carpetas para los diferentes escenarios/anos en el futuro donde se har?n los modelos
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1")#COLOCAR AQUI LA RUTA DEL DIRECTORIO DE DONDE Guardar las capas
dir.create("MH_CCSM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("MH_MIROC")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("MH_MPI_ESM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LGM_CCSM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LGM_MIROC")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LGM_MPI_ESM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LIG_CCSM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("M_total_presente")###comando para crear una carpeta con el nombre indicado entre comillas

##PART II: Cut the variable to the "M" extension
EME2 <- readOGR("D:/proyectos_sigs/grasshoppers_project/m_finales/S_zapotecum.shp")
plot(EME2)##Plot the individual "M" to species

EME <- readOGR("D:/proyectos_sigs/grasshoppers_project/m_finales/S_total.shp")
plot(EME)##Plot the combined "M" to species group.



##Bioclimate information to PRESENT
#PRESENT_M_species
#setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/presente/")
setwd("D:/proyectos_sigs/salomon_diversiﬁcation of grasshoppers/modelos_pasado/kuenm_models_spp/Sphenarium_rugosum/G_variables/Set_1/M_total_presente/")
pca_path <- list.files(".",pattern = "*.asc",full.names = T)###crea el stack de los COMPONENTES 
capas_presente<- stack(pca_path)

##PRESENT_M_species
CROPPED_EME2 <- crop(capas_presente, extent(EME2))
MASKED_EME2 <- mask(CROPPED_EME2, EME2)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME2), function(x){
  writeRaster(MASKED_EME2[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/M_variables/Set_1/", x,".asc"),overwrite=TRUE)})


##PRESENT_M_total
CROPPED_EME <- crop(capas_presente, extent(EME))
MASKED_EME <- mask(CROPPED_EME, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME), function(x){
  writeRaster(MASKED_EME[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/M_total_presente/", x,".asc"),overwrite=TRUE)})


#HOLOCENO: MH_CCSM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/MH_CCSM/")
pca_path2 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_MH_CCSM<- stack(pca_path2)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_MH_CCSM <- crop(capas_MH_CCSM, extent(EME))
MASKED_EME_MH_CCSM <- mask(CROPPED_EME_MH_CCSM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_MH_CCSM), function(x){
  writeRaster(MASKED_EME_MH_CCSM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/MH_CCSM/", x,".asc"),overwrite=TRUE)})

#HOLOCENO: MH_MIROC
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/MH_MIROC/")
pca_path3 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_MH_MIROC<- stack(pca_path3)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_MH_MIROC <- crop(capas_MH_MIROC, extent(EME))
MASKED_EME_MH_MIROC <- mask(CROPPED_EME_MH_MIROC, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_MH_MIROC), function(x){
  writeRaster(MASKED_EME_MH_MIROC[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/MH_MIROC/", x,".asc"),overwrite=TRUE)})

#HOLOCENO: MH_MPI_ESM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/MH_MPI_ESM/")
pca_path4 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_MH_MPI_ESM<- stack(pca_path4)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_MH_MPI_ESM <- crop(capas_MH_MPI_ESM, extent(EME))
MASKED_EME_MH_MPI_ESM <- mask(CROPPED_EME_MH_MPI_ESM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_MH_MPI_ESM), function(x){
  writeRaster(MASKED_EME_MH_MPI_ESM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/MH_MPI_ESM/", x,".asc"),overwrite=TRUE)})


#Last Maximun Glacial: LGM_CCSM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/LGM_CCSM/")
pca_path5 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_LGM_CCSM<- stack(pca_path5)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LGM_CCSM <- crop(capas_LGM_CCSM, extent(EME))
MASKED_EME_LGM_CCSM <- mask(CROPPED_EME_LGM_CCSM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LGM_CCSM), function(x){
  writeRaster(MASKED_EME_LGM_CCSM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LGM_CCSM/", x,".asc"),overwrite=TRUE)})


#Last Maximun Glacial: LGM_MIROC
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/LGM_MIROC/")
pca_path6 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_LGM_MIROC<- stack(pca_path6)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LGM_MIROC <- crop(capas_LGM_MIROC, extent(EME))
MASKED_EME_LGM_MIROC <- mask(CROPPED_EME_LGM_MIROC, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LGM_MIROC), function(x){
  writeRaster(MASKED_EME_LGM_MIROC[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LGM_MIROC/", x,".asc"),overwrite=TRUE)})


#Last Maximun Glacial: LGM_MPI_ESM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/LGM_MPI_ESM/")
pca_path7 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_LGM_MPI_ESM<- stack(pca_path7)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LGM_MPI_ESM <- crop(capas_LGM_MPI_ESM, extent(EME))
MASKED_EME_LGM_MPI_ESM <- mask(CROPPED_EME_LGM_MPI_ESM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LGM_MPI_ESM), function(x){
  writeRaster(MASKED_EME_LGM_MPI_ESM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LGM_MPI_ESM/", x,".asc"),overwrite=TRUE)})


#Last Inter-Glacial: LIG_CCSM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/lig_bios_final/")
pca_path8 <- list.files(".",pattern = "*.asc",full.names = T)###crea el stack de los COMPONENTES 
capas_LIG_CCSM<- stack(pca_path8)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LIG_CCSM <- crop(capas_LIG_CCSM, extent(EME))
MASKED_EME_LIG_CCSM <- mask(CROPPED_EME_LIG_CCSM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LIG_CCSM), function(x){
  writeRaster(MASKED_EME_LIG_CCSM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LIG_CCSM/", x,".asc"),overwrite=TRUE)})



rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo
###############################################
###############################################
##Correr los modelos en la plataforma ENMval###
###############################################
###############################################
#establecemos la carpeta de trabajo en el directorio original
options(java.parameters = "-Xmx8000m")
library(ENMeval)
library(raster)
library(dplyr)
library(rJava)
library(kuenm)

# Set a random seed in order to be able to reproduce this analysis.
set.seed(48)

#Llamar el archivo CSV final que tiene los datos definitivos de la especie
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/")
datos <- read.csv2("Sphenarium_zapotecum.csv", sep = ",", header = TRUE)## leer el archivo de puntos de inter?s
names(datos)##vemos como se llaman las columnas de la tabla de datos
datos$lat <- as.numeric(datos$lat)###volver las variables en formato n?merico
datos$lon <- as.numeric(datos$lon)###volver las variables en formato n?merico

occs <- (datos[, 2:3])
occs <- occs[!duplicated(occs),]# avoid pseudoreplication.


##indicar donde estan las carpetas del presente de la especie a modelar
envs <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/M_variables/Set_1", full.names = TRUE))

##crear los datos de background para poder generar luego los modelos
occs.cells <- raster::extract(envs[[1]], occs[,1:2], cellnumbers = TRUE)
occs.cellDups <- duplicated(occs.cells[,1])
occs <- occs[!occs.cellDups,]

occs.sf <- sf::st_as_sf(occs, coords = c("lon","lat"), crs = raster::crs(envs))
occs.sf <- sf::st_as_sf(occs, coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#Hacer lo mismo con las capas climaticas
crs(envs) <- raster::crs(occs.sf)


occs.buf <- sf::st_buffer(occs.sf, dist = 120000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))

# To add sf objects to a plot, use add = TRUE
plot(envs[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)


# Crop environmental rasters to match the study extent
envs.bg <- raster::crop(envs, occs.buf)
envs.bg <- raster::mask(envs.bg, occs.buf)

plot(envs.bg[[1]], main = names(envs)[1])


# other rasters have data, ENMeval internally converts these cells to NA.
bg <- dismo::randomPoints(envs.bg[[1]], n = 20000) %>% as.data.frame()
colnames(bg) <- colnames(occs)

# Notice how we have pretty good coverage (every cell).
plot(envs.bg[[1]])
points(bg, pch = 20, cex = 0.2)


#Correr el maxent bajo el supuesto de Jackknife.
##Activar la funci?n Partial ROC para el calculo de la metrica
proc <- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(proc_auc_ratio = proc$pROC_summary[1], 
                    proc_pval = proc$pROC_summary[2], row.names = NULL)
  return(out)
}


##model1
e.mx.l <- ENMevaluate(occs = occs, envs = c(envs.bg), bg = bg, 
                      algorithm = 'maxent.jar', 
                      partitions = 'jackknife',
                      tune.args = list(fc = c("L","LQ","LQH","H","Q", "T", "LPT"), rm = 1:8), user.eval = proc)

e.mx.l@results
e.mx.l@results.partitions

#Seleccionar al mejor modelo de acuerdo al AICc
res <- eval.results(e.mx.l)
res.partitions <- eval.results.partitions(e.mx.l)
opt.aicc <-  res %>%
  dplyr::filter(!is.na(AICc)) %>% 
  dplyr::filter(!is.na(proc_auc_ratio.avg)) %>% 
  dplyr::filter(AICc == min(AICc)) %>%
  dplyr::filter(proc_auc_ratio.avg == max(proc_auc_ratio.avg )) %>%
  slice(1)

output_table <- data.frame(Especie = "Sphenarium_zapotecum",
                           RM = opt.aicc$rm,
                           FC = opt.aicc$fc,
                           AICc = opt.aicc$AICc,
                           ROC_Partial = opt.aicc$proc_auc_ratio.avg)

###
Variable_importance <- opt.aicc$tune.args %>% purrr::map(function(x){
  
  e.mx.l@variable.importance[[x]]   %>% mutate(model = x) %>% as_tibble()
}) %>% purrr::reduce(bind_rows)

percent.contribution_df <- Variable_importance %>%
  select(-percent.contribution) %>% 
  tidyr::spread(variable, permutation.importance)
percent.contribution_df

mod.seq <- eval.models(e.mx.l)[[opt.aicc$tune.args]]
presente <- eval.predictions(e.mx.l)[[opt.aicc$tune.args]]
plot(presente)

##proyectar los resultados del maxent en los mapas para crear las predicciones (presente y futuro)
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum")#COLOCAR AQUI LA RUTA DEL DIRECTORIO DE DONDE Guardar las capas
dir.create("Model_results")##
writeRaster(x = presente, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_m_spp.asc", overwrite=TRUE)#guarda el mapa del presente

##Guardar los archivos CSV que tienen la informaci?n de los mejores modelos.
write.csv(x = opt.aicc, file = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/best_model_set1.csv")
write.csv(x = res.partitions, file = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/output_results_partitions_set1.csv")
write.csv(x = Variable_importance, file = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/variables_importance_set1.csv")


##indicar donde estan las carpetas de las capas climaticas del pasado para la especie a modelar
shp_presente <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/M_total_presente/", full.names = TRUE))##ano tipico 2
shp_holoceno1 <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/MH_CCSM", full.names = TRUE))##ano tipico 3
shp_holoceno2 <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/MH_MIROC", full.names = TRUE))##ano tipico 4
shp_holoceno3 <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/MH_MPI_ESM", full.names = TRUE))##ano tipico 5
shp_LGM1 <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LGM_CCSM", full.names = TRUE))##ano tipico 1
shp_LGM2 <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LGM_MIROC", full.names = TRUE))##ano tipico 2
shp_LGM3 <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LGM_MPI_ESM", full.names = TRUE))##ano tipico 3
shp_LIG1 <- raster::stack(list.files("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/G_variables/Set_1/LIG_CCSM", full.names = TRUE))##ano tipico 4


##proyectar los resultados del maxent en los mapas para crear las predicciones (presente y futuro)
FY1 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_presente)##mapa del FY1
FY2 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_holoceno1)##mapa del FY2
FY3 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_holoceno2)##mapa del FY3
FY4 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_holoceno3)##mapa del FY4
FY5 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_LGM1)##mapa del FY5
FY6 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_LGM2)##mapa del FY5
FY7 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_LGM3)##mapa del FY5
FY8 <- predict(e.mx.l@models[[opt.aicc$tune.args]], shp_LIG1)##mapa del FY5


setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/")#COLOCAR AQUI LA RUTA DEL DIRECTORIO DE DONDE Guardar las capas

writeRaster(x = FY1, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_m_total.asc", overwrite=TRUE)#guarda el mapa del presente
writeRaster(x = FY2, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_holoceno1.asc", overwrite=TRUE)#guarda el mapa del presente
writeRaster(x = FY3, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_holoceno2.asc", overwrite=TRUE)#guarda el mapa del presente
writeRaster(x = FY4, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_holoceno3.asc", overwrite=TRUE)#guarda el mapa del presente
writeRaster(x = FY5, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LGM1.asc", overwrite=TRUE)#guarda el mapa del presente
writeRaster(x = FY6, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LGM2.asc", overwrite=TRUE)#guarda el mapa del presente
writeRaster(x = FY7, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LGM3.asc", overwrite=TRUE)#guarda el mapa del presente
writeRaster(x = FY8, filename = "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LIG1.asc", overwrite=TRUE)#guarda el mapa del presente


########################################################
##########ANALISIS ESPACIALES DE MODELOS################
########################################################
rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo

########################################################################
######Primera PARTE: Calcular el mapa binario del presente #############
########################################################################
#1.llamar el modelo promedio que se obtivo para el presente
presente <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_m_spp.asc")
plot(presente)

#2. Llamar al archivo .csv de los puntos utilizados para el training del modelo final.
data_spp <- read.csv2("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Sphenarium_zapotecum.csv", sep = ",", header = TRUE)
data_spp$lat <- as.numeric(as.character(data_spp$lat))###volver las variables en formato n?merico
data_spp$lon <- as.numeric(as.character(data_spp$lon))###volver las variables en formato n?merico

#1. Convertir los puntos en datos de presencia para un SIG
points_occ_fin <- SpatialPointsDataFrame(data_spp[,2:3],data_spp)#convertir a un archivo shp de puntos 

#3. Extraer los valores de idoneidad que el modelo predice para cada uno de los puntos de presencia (trainning)
presencia_model <- data.frame(raster::extract(presente, points_occ_fin [,2:3]))###extrae los valorespara cada uno de los puntos
presencia_model2 <- na.omit(presencia_model)## omite mis datos de presencia sin valores climaticos

#4. calcular el valor del 10% de los datos de presencia de la especie
umbral_models <- quantile(presencia_model2$raster..extract.presente..points_occ_fin...2.3.., prob= 0.09)##me indica el valor del umbral

#5. Transformar la capa del presente a un mapa binario
presente_bin <- presente >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia

plot(presente_bin)###me grafica el mapa
plot(points_occ_fin, add= T, pch=19, col= "red")###grafica los puntos de entrenamiento

###Salvar el mapa de la especie
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/presente_M/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(presente_bin, filename= "Sphenarium_zapotecum.asc", overwrite=T, suffix='names')


######################
##PRESENTE todo#######
######################
#1.llamar el modelo promedio que se obtivo para el presente
presente_M_total <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_m_total.asc")
plot(presente_M_total)
#5. Transformar la capa del presente a un mapa binario
presente_M_total_bin1 <- presente_M_total >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
plot(presente_M_total_bin1)


###Salvar el mapa de la especie
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/presente_M_total/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(presente_M_total_bin1, filename= "Sphenarium_zapotecum.asc", overwrite=T, suffix='names')


#################
##HOLOCENO#######
#################
#1.llamar el modelo promedio que se obtivo para el presente
Holoceno1 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_holoceno1.asc")
Holoceno2 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_holoceno2.asc")
Holoceno3 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_holoceno3.asc")

#5. Transformar la capa del presente a un mapa binario
holoceno_bin1 <- Holoceno1 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
holoceno_bin2 <- Holoceno2 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
holoceno_bin3 <- Holoceno3 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia

holoceno_fin <- holoceno_bin1 + holoceno_bin2 + holoceno_bin3
plot(holoceno_fin)
holoceno_fin_bin <- holoceno_fin >= 2
plot(holoceno_fin_bin)


###Salvar el mapa de la especie
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/holoceno/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(holoceno_fin_bin, filename= "Sphenarium_zapotecum.asc", overwrite=T, suffix='names')


#################
##LGM#######
#################
#1.llamar el modelo promedio que se obtivo para el presente
LGM1 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LGM1.asc")
LGM2 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LGM2.asc")
LGM3 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LGM3.asc")

#5. Transformar la capa del presente a un mapa binario
LGM_bin1 <- LGM1 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
LGM_bin2 <- LGM2 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
LGM_bin3 <- LGM3 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia

LGM_fin <- LGM_bin1 + LGM_bin2 + LGM_bin3
plot(LGM_fin)
LGM_fin_bin <- LGM_fin >= 2
plot(LGM_fin_bin)


###Salvar el mapa de la especie
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/LGM/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(LGM_fin_bin, filename= "Sphenarium_zapotecum.asc", overwrite=T, suffix='names')


######################
##LIG todo#######
######################
#1.llamar el modelo promedio que se obtivo para el presente
LIG <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_zapotecum/Model_results/Sphenarium_zapotecum_LIG1.asc")

#5. Transformar la capa del presente a un mapa binario
LIG_bin1 <- LIG >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
plot(LIG_bin1)


###Salvar el mapa de la especie
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/LIG/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(LIG_bin1, filename= "Sphenarium_zapotecum.asc", overwrite=T, suffix='names')


