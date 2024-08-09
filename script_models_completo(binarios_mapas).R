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
library(ecospat)
library(dplyr)
library(reshape)
library(kuenm)

rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo
##################################################
##CUT THE BIOCLIMATE VARIABLES INTO "M" EXTENT####
##################################################
#PARTE I: Create the folders and directories
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/")#
dir.create("Sphenarium_rugosum")###comando para crear una carpeta con el nombre indicado entre comillas

setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/")
dir.create("G_variables")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("M_variables")###comando para crear una carpeta con el nombre indicado entre comillas

#3.Crear las carpetas para los dos "Set" diferentes de variables del presente (Set1)
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/M_variables/")#
dir.create("Set_1")###comando para crear una carpeta con el nombre indicado entre comillas

#4.Crear las carpetas para los dos "Set" diferentes de variables del futuro (Set1)
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/")#
dir.create("Set_1")###comando para crear una carpeta con el nombre indicado entre comillas

#5.Crear las carpetas para los diferentes escenarios/anos en el futuro donde se har?n los modelos
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1")#COLOCAR AQUI LA RUTA DEL DIRECTORIO DE DONDE Guardar las capas
dir.create("MH_CCSM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("MH_MIROC")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("MH_MPI_ESM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LGM_CCSM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LGM_MIROC")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LGM_MPI_ESM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("LIG_CCSM")###comando para crear una carpeta con el nombre indicado entre comillas
dir.create("M_total_presente")###comando para crear una carpeta con el nombre indicado entre comillas

##PART II: Cut the variable to the "M" extension
EME2 <- readOGR("D:/proyectos_sigs/grasshoppers_project/m_finales/S_rugosum.shp")
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
  writeRaster(MASKED_EME2[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/M_variables/Set_1/", x,".asc"),overwrite=TRUE)})


##PRESENT_M_total
CROPPED_EME <- crop(capas_presente, extent(EME))
MASKED_EME <- mask(CROPPED_EME, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME), function(x){
  writeRaster(MASKED_EME[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/M_total_presente/", x,".asc"),overwrite=TRUE)})


#HOLOCENO: MH_CCSM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/MH_CCSM/")
pca_path2 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_MH_CCSM<- stack(pca_path2)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_MH_CCSM <- crop(capas_MH_CCSM, extent(EME))
MASKED_EME_MH_CCSM <- mask(CROPPED_EME_MH_CCSM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_MH_CCSM), function(x){
  writeRaster(MASKED_EME_MH_CCSM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/MH_CCSM/", x,".asc"),overwrite=TRUE)})

#HOLOCENO: MH_MIROC
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/MH_MIROC/")
pca_path3 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_MH_MIROC<- stack(pca_path3)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_MH_MIROC <- crop(capas_MH_MIROC, extent(EME))
MASKED_EME_MH_MIROC <- mask(CROPPED_EME_MH_MIROC, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_MH_MIROC), function(x){
  writeRaster(MASKED_EME_MH_MIROC[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/MH_MIROC/", x,".asc"),overwrite=TRUE)})

#HOLOCENO: MH_MPI_ESM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/MH_MPI_ESM/")
pca_path4 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_MH_MPI_ESM<- stack(pca_path4)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_MH_MPI_ESM <- crop(capas_MH_MPI_ESM, extent(EME))
MASKED_EME_MH_MPI_ESM <- mask(CROPPED_EME_MH_MPI_ESM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_MH_MPI_ESM), function(x){
  writeRaster(MASKED_EME_MH_MPI_ESM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/MH_MPI_ESM/", x,".asc"),overwrite=TRUE)})


#Last Maximun Glacial: LGM_CCSM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/LGM_CCSM/")
pca_path5 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_LGM_CCSM<- stack(pca_path5)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LGM_CCSM <- crop(capas_LGM_CCSM, extent(EME))
MASKED_EME_LGM_CCSM <- mask(CROPPED_EME_LGM_CCSM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LGM_CCSM), function(x){
  writeRaster(MASKED_EME_LGM_CCSM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/LGM_CCSM/", x,".asc"),overwrite=TRUE)})


#Last Maximun Glacial: LGM_MIROC
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/LGM_MIROC/")
pca_path6 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_LGM_MIROC<- stack(pca_path6)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LGM_MIROC <- crop(capas_LGM_MIROC, extent(EME))
MASKED_EME_LGM_MIROC <- mask(CROPPED_EME_LGM_MIROC, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LGM_MIROC), function(x){
  writeRaster(MASKED_EME_LGM_MIROC[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/LGM_MIROC/", x,".asc"),overwrite=TRUE)})


#Last Maximun Glacial: LGM_MPI_ESM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/LGM_MPI_ESM/")
pca_path7 <- list.files(".",pattern = "*.tif",full.names = T)###crea el stack de los COMPONENTES 
capas_LGM_MPI_ESM<- stack(pca_path7)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LGM_MPI_ESM <- crop(capas_LGM_MPI_ESM, extent(EME))
MASKED_EME_LGM_MPI_ESM <- mask(CROPPED_EME_LGM_MPI_ESM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LGM_MPI_ESM), function(x){
  writeRaster(MASKED_EME_LGM_MPI_ESM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/LGM_MPI_ESM/", x,".asc"),overwrite=TRUE)})


#Last Inter-Glacial: LIG_CCSM
setwd("D:/proyectos_sigs/grasshoppers_project/capas_climaticas/lig_bios_final/")
pca_path8 <- list.files(".",pattern = "*.asc",full.names = T)###crea el stack de los COMPONENTES 
capas_LIG_CCSM<- stack(pca_path8)

#Recortar las capas a la forma de la "M" de la especie de inter?s.
CROPPED_EME_LIG_CCSM <- crop(capas_LIG_CCSM, extent(EME))
MASKED_EME_LIG_CCSM <- mask(CROPPED_EME_LIG_CCSM, EME)

#Funcion para mandar a guardar de forma automatica las capas recortadas:
lapply(names(MASKED_EME_LIG_CCSM), function(x){
  writeRaster(MASKED_EME_LIG_CCSM[[x]], paste0("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/G_variables/Set_1/LIG_CCSM/", x,".asc"),overwrite=TRUE)})


###############################################
##Correr los modelos en la plataforma KUENM####
###############################################
###############################################
rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo para evitar confusiones en archivos

###################################################################################################
#Parte 1: Dividir en los datos en entrenamiento y validacion para hacer el protocolo de modelado###
###################################################################################################
#Llamar el archivo CSV final que tiene los datos definitivos de la especie
data<- read.csv2("D:/proyectos_sigs/grasshoppers_project/Sphenarium_rugosum.csv", sep = ",", header = TRUE)## leer el archivo de puntos de inter?s
data$lat <- as.numeric(data$lat)###volver las variables en formato n?merico
data$lon <- as.numeric(data$lon)###volver las variables en formato n?merico

todos <- unique(data)###seleccionar el archivo base donde estan todos los datos
todos$sp <- paste(todos[,2], todos[,3])###seleccionar la columna lon y lat de los datos
train <- todos[sample(nrow(todos), round((length(todos[,1])/4 *3))), ]###indicar como se hará la seleccion de los datos (en este caso es seleccionar 4 de las 5 partes = 80%)
test_ind <- todos[!todos[,4] %in% train[,4], ]### coloca una nueva columna donde se indica la seleccion realizada
#los registros han sido dividos en tres objetos, todos (BD original), train (80% del total de la BD original) y test que corresponde al 20% de datos individuales de evaluación (20% del total de la BD original)

#guardamos nuestros subconjuntos
#Datos de testing Independientes para evaluar el modelo final
write.csv(test_ind[,1:3], "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Sphenarium_rugosum_ind.csv", row.names = FALSE)###nombre con el que guardará el archivo .csv que contiene el 20% de datos independientes de evaluacion

#Datos de training juntos (80%) para construir el modelo final
write.csv(train[,1:3], "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Sphenarium_rugosum_joint.csv", row.names = FALSE)###nombre con el que guardará el archivo .csv que contiene el 80% de datos originales


#5. Dividir el archivo de los datos de trainning en 2 archivos (training [calibration]vs. testing [calibration]) para la calibracion del modelo
rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo
data1<- read.csv2("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Sphenarium_rugosum_joint.csv", sep=",", header=TRUE) ###seleccionar el archivo que ocntine los datos de presencia de la especie a modelar
data1$lat <- as.numeric(as.character(data1$lat))##convertir la informaci?n en n?mero
data1$lon <- as.numeric(as.character(data1$lon))##convertir la informaci?n en n?mero

todos <- unique(data1)###seleccionar el archivo base donde estan todos los datos
todos$sp <- paste(todos[,2], todos[,3])###seleccionar la columna lon y lat de los datos
train <- todos[sample(nrow(todos), round((length(todos[,1])/4 *3))), ]###indicar como se hará la seleccion de los datos (en este caso es seleccionar 4 de las 5 partes = 80%)
test_test <- todos[!todos[,4] %in% train[,4], ]### coloca una nueva columna donde se indica la seleccion realizada
#los registros han sido dividos en tres objetos

#guardamos nuestros subconjuntos
write.csv(train[,1:3], "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Sphenarium_rugosum_train.csv", row.names = FALSE)###nombre con el que guardara el archivo .csv que contiene el 80% de datos para calibrar los modelos

write.csv(test_test[,1:3], "D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Sphenarium_rugosum_test.csv", row.names = FALSE)###nombre con el que guardara el archivo .csv que contiene el 20% de los datos para la calibraci?n del modelo



#######################################################
######### Parte 2: Correr el modelo ###################
#######################################################
rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo
############################################
######### ETAPA 1 #########################
setwd("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/")### directorio donde estan todos los archivos 

############################################
######### The model creation ###############
############################################
occ_joint <- "Sphenarium_rugosum_joint.csv"##archivo que tiene los datos completos
occ_tra <- "Sphenarium_rugosum_train.csv"##archivo que tiene el 80% de los datos para entrenar
M_var_dir <- "M_variables"##carpeta donde estan las capas del presente
batch_cal <- "Candidate_models"##carpeta donde se guardaran los resultados de los 620 modelos candidatos
out_dir <- "Candidate_Models"##carpeta donde se guardaran los resultados de los 620 modelos candidatos
reg_mult <- c(seq(0.2, 2, 0.2), seq(2, 6, 0.5),8, 10)##indicar todas las combinaciones de RM que se usaran
f_clas <- "all"###indicar cuales combinaciones de "Feature" se utilizara
args <- NULL 

#Correr el programa para que haga 620 modelos simultaneos que leugo se usaran para definir al mejor.
maxent_path <- "C:/Users/ASUS/Documents/R/win-library/4.1/dismo/java" ##esto debe localizarlo en su computadora
wait <- FALSE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)
###Se abrir? una ventana negra con letras blancas que ir? indicando como va el avance para la elaboraci?n de los 620 modelos... esa ventana se cierra sola... cuando se cierre significa que ya termino y puede seguir ejecutando las otras partes del script.ESTO PUEDE TARDAR VARIOS MINUTOS E INCLUSO HORAS (DEPENDE DE LA CAPACIDAD DEL COMPUTADOR)...

############################################
######### The model evaluation #############
############################################
occ_test <- "Sphenarium_rugosum_test.csv"##indicar donde esta el archivo para evaluar a los 620 modelos
out_eval <- "Calibration_results"##indicar donde guardar los resultados a obtener
threshold <- 5##cual es el criterio de seleccion en referencia al valor de error de omision
rand_percent <- 50##porcentaje a utilizar de los datos de forma aleatoria
iterations <- 500##numero de iteracciones
kept <- TRUE
selection <- "OR_AICc"##criterio de seleccion de modelos: omission error and Akaike criteria
paral_proc <- FALSE ##no correr en paralelo


cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)

###Aparecera una ventana en blanco con una barra de progreso.. esa barra se ir? colocando de color azul a medida que el programa avance en el an?lisis... al terminar (llegar al 100%) se cerrar? y entonces podr? avanzar en el script. ESTO IGUAL PUEDE TARDAR VARIOS MINUTOS Y/O HORAS, recuerde que esta evaluando 620 modelos (PERO NO es normal que no avance nada en m?s de 2 horas... si eso ocurre hay alg?n problema y debe revisarse).


############################################
######### The FINAL model creation ######### 
############################################
batch_fin <- "Final_models"#donde guardar el(los) modelos finales que se seleccionaron
mod_dir <- "Final_Models"#donde guardar el(los) modelos finales que se seleccionaron
rep_n <- 10#indicar cuantas replicas se haran
rep_type <- "Bootstrap"##cual es el criterio para hacer replicas
jackknife <- FALSE
out_format <- "cloglog"
project <- TRUE##indicar que queremos que haga modelos a futuro
G_var_dir <- "G_variables"##indicar donde estan las capas del futuro
ext_type <- "ext_clam"#indicar que queremos que permita extrapolacion y docampling
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL 
kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)
###Se abrira una ventana negra con letras blancas que ira indicando como va el avance para la elaboracion de (los) modelos finales... este sera el modelo utilizando los mejores parametros de configuracion esa ventana se cierra sola... cuando se cierre significa que ya termino y puede seguir ejecutando las otras partes del script.
#esto puede tardar un poco (quizas un par de horas)... depende de la PC, la cantidad de datos y el tamano del area


############################################
####### The model evaluation final #########
############################################
occ_ind <- "Sphenarium_rugosum_ind.csv"##Elegir el archivo de 20% de los datos que nunca se usaron para evaluar el modelo
replicates <- TRUE##indicar que hay replicas
threshold <- 10
out_feval <- "Final_Models_evaluation"##nombre del archivo donde se guardaran los resultados

fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)

best <- read.csv("Calibration_results/selected_models.csv")
knitr::kable(best, caption = "Models selected based on significance, omission rates, and AICc, in that order.")
##al terminar de leer esta parte le saldra una tabla final con resultados... esos resultados deben ser copiados en el archivo excel, llenando cada una de las cosas que se solitica.. Si en la tabla le sale m?s de una opci?n (como es este ejemplo) ustede debe solo copiar los datos del primer modelo... Debe copiar de ac? los valores de RM, feature y AICc


best <- read.csv("Final_Models_evaluation/fm_evaluation_results.csv")
knitr::kable(best, caption = "Models selected based on significance, omission rates, and AICc, in that order.")


########################################################
##########ANALISIS ESPACIALES DE MODELOS################
########################################################
rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo

########################################################################
######Primera PARTE: Calcular el mapa binario del presente #############
########################################################################
#1.llamar el modelo promedio que se obtivo para el presente
presente <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_median.asc")
plot(presente)

#2. Llamar al archivo .csv de los puntos utilizados para el training del modelo final.
data_spp <- read.csv2("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Sphenarium_rugosum_joint.csv", sep = ",", header = TRUE)
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

#6. Llamar al archivo .csv de los puntos utilizados para el testing del modelo final (es el archivo que se llama "IND").
data_testing_spp <- read.csv2("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Sphenarium_rugosum_ind.csv", sep = ",", header = TRUE)
data_testing_spp$lat <- as.numeric(as.character(data_testing_spp$lat))###volver las variables en formato n?merico
data_testing_spp$lon <- as.numeric(as.character(data_testing_spp$lon))###volver las variables en formato n?merico

#7. Convertir los puntos en datos de presencia para un SIG
points_test_fin <- SpatialPointsDataFrame(data_testing_spp[,2:3],data_testing_spp)#convertir a un archivo shp de puntos 

#8. Extraer los valores de idoneidad que el modelo predice para cada uno de los puntos de presencia (Testing)
presencia_model_test <- data.frame(raster::extract(presente_bin, points_test_fin [,2:3]))###extrae los valores para saber cuantos datos del TESTING son predichos correctamente y cuantos no.
presencia_model_test2 <- na.omit(presencia_model_test)## omite mis datos de presencia sin valores climaticos
sum(presencia_model_test2$raster..extract.presente_bin..points_test_fin...2.3..)###me indica en n?mero cuantos datos si fueron bien predichos


###Salvar el mapa de la especie
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/presente_M/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(presente_bin, filename= "Sphenarium_rugosum.asc", overwrite=T, suffix='names')


presente_mtotal <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_M_total_presente_median.asc")##modelo del CCSM
presente_mtotal_bin <- presente_mtotal >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
plot(presente_mtotal_bin)
#4.Guardar el mapa de Presente_total 
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/presente_M_total/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(presente_mtotal_bin, filename= "Sphenarium_rugosum.asc", overwrite=T, suffix='names')


########################################################################
######Segunda PARTE: Calcular el mapa binario del HOLOCENO##############
########################################################################
#1.llamar las capas (modelos/replicas) que se hicieron del futuro 2040 (uno por cada laboratorio)
Holoceno_1 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_MH_CCSM_median.asc")##modelo del CCSM
Holoceno_2 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_MH_MIROC_median.asc")##modelo del MIROC
Holoceno_3 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_MH_MPI_ESM_median.asc")##modelo del MPI_ESM

#2.Transformar los mapas a mapas binarios
Holoceno_1_bin <- Holoceno_1 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
Holoceno_2_bin <- Holoceno_2 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
Holoceno_3_bin <- Holoceno_3 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia

#3. Sumar los 5 modelos para obtener el ?rea de coincidencia 
map_holoceno_total <- Holoceno_1_bin + Holoceno_2_bin + Holoceno_3_bin 
plot(map_holoceno_total)
map_holoceno_fin <- map_holoceno_total >= 2 ##selecciona como "presencia" aquellas ?reas donde al menos 3 modelos coinciden.

#4.Guardar el mapa
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/holoceno/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(map_holoceno_fin, filename= "Sphenarium_rugosum.asc", overwrite=T, suffix='names')


###################################################################
######Tercera PARTE: Calcular el mapa binario del LGM##############
###################################################################
#1.llamar las capas (modelos/replicas) que se hicieron del futuro 2040 (uno por cada laboratorio)
LGM_1 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_LGM_CCSM_median.asc")##modelo del CCSM
LGM_2 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_LGM_MIROC_median.asc")##modelo del MIROC
LGM_3 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_LGM_MPI_ESM_median.asc")##modelo del MPI_ESM

#2.Transformar los mapas a mapas binarios
LGM_1_bin <- LGM_1 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
LGM_2_bin <- LGM_2 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
LGM_3_bin <- LGM_3 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia

#3. Sumar los 5 modelos para obtener el ?rea de coincidencia 
map_LGM_total <- LGM_1_bin + LGM_2_bin + LGM_3_bin 
plot(map_LGM_total)
map_LGM_fin <- map_LGM_total >= 2 ##selecciona como "presencia" aquellas ?reas donde al menos 3 modelos coinciden.

#4.Guardar el mapa
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/LGM/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(map_LGM_fin, filename= "Sphenarium_rugosum.asc", overwrite=T, suffix='names')


###################################################################
######Cuarta PARTE: Calcular el mapa binario del LIG##############
###################################################################
#1.llamar las capas (modelos/replicas) que se hicieron del futuro 2040 (uno por cada laboratorio)
LIG_1 <- raster("D:/proyectos_sigs/grasshoppers_project/modelos_pasado/Sphenarium_rugosum/Final_Models/M_1.2_F_lqp_Set_1_EC/Sphenarium_rugosum_LIG_CCSM_median.asc")##modelo del CCSM


#2.Transformar los mapas a mapas binarios
LIG_1_bin <- LIG_1 >= umbral_models #convierte el mapa de idoneidad en presencia/ausencia
plot(LIG_1_bin)

#4.Guardar el mapa 
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/LIG/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(LIG_1_bin, filename= "Sphenarium_rugosum.asc", overwrite=T, suffix='names')

################################################
######Quinta PARTE: Estables areas##############
################################################
#3. Sumar los 5 modelos para obtener el ?rea de coincidencia 
mapa_estable_areas <- presente_bin + map_holoceno_fin + map_LGM_fin + LIG_1_bin 
plot(mapa_estable_areas)
mapa_estable_areas_fin <- mapa_estable_areas >= 4 ##selecciona como "presencia" aquellas ?reas donde al menos 3 modelos coinciden.
plot(mapa_estable_areas_fin)

#4.Guardar el mapa
setwd("D:/proyectos_sigs/grasshoppers_project/mapas_binarios/estables_areas/")#RUTA DEL DIRECTORIO DE DONDE Guardar las capas
writeRaster(mapa_estable_areas_fin, filename= "Sphenarium_rugosum.asc", overwrite=T, suffix='names')


raster::freq(presente_bin)#en este caso hay 37596 pixeles con valor 1 (presencia) este valor se coloca en el archivo excel 
raster::freq(map_holoceno_fin)#en este caso hay 19737 pixeles con valor 1 (presencia) este valor se coloca en el archivo excel 
raster::freq(map_LGM_fin)#en este caso hay 19737 pixeles con valor 1 (presencia) este valor se coloca en el archivo excel 
raster::freq(LIG_1_bin)#en este caso hay 19737 pixeles con valor 1 (presencia) este valor se coloca en el archivo excel 
raster::freq(mapa_estable_areas_fin)#en este caso hay 19737 pixeles con valor 1 (presencia) este valor se coloca en el archivo excel 
##FIN!
