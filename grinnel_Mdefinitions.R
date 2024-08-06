## Dependences
library(terra)
library(grinnell)

setwd("D:/diversiﬁcation_grasshoppers/M_analyses/")##indicate the main rute to work

#Species ocurrence data
sp_files <- read.csv("./Sphenarium_variabile/occurrences.csv", header = T)

#Bioclimate variable data
variables <- terra::rast(list.files("D:/diversiﬁcation_grasshoppers/capas_climaticas/presente/",
                                      pattern = ".asc$", full.names = TRUE))
plot(variables$bio_1)##draw and review the first climate variable

# removing records that are outside variables (only to confirm the quality of information into data)
records <- sp_files[!is.na(raster::extract(variables[[1]], sp_files[, 2:3])), ]
records2 <- na.omit(records)

#Running in current scenario
m <- M_simulationR(data = records2, current_variables = variables,
                   dispersal_kernel = "normal", kernel_spread = 5,
                   max_dispersers = 10, replicates = 16, dispersal_events = 50,
                   output_directory = "S_variabile_m")



