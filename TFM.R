################################################

# Script TFM - UEM

#Alumno: Sandra Roldan 
###############################################

############################
### INSTALACION PAQUETES ###
############################

if (!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager")}
if (!requireNamespace("edgeR", quietly=TRUE)){BiocManager::install("edgeR")}
if (!requireNamespace("tidyverse", quietly=TRUE)){install.packages("tidyverse")}
BiocManager::install("hipathia")


###############################
### PREPARACION ENVIRONMENT ###
###############################

library(edgeR)
library(BiocManager)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(genefu)
library(tidyverse)
library(hipathia)


################################################################
####  Importacion de archivos y filtracion/cruzado de datos ####
################################################################

#Matriz de datos de expresion:
genexp <- read.delim(file="./Raw_Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2)



#Dimensiones de las matrices
dim(genexp)
#genexp1<-genexp[,-1]   #elimina la columna name
genexp1<-genexp[,-2]   #elimina la columna description


#Comprobar si hay duplicados
#sum(duplicated(genexp1$Description))
sum(duplicated(genexp1$Name))


#Cambiamos los nombres de fila por los nombres de muestra
#variable name: nomenclatura de genes "Ensembl"
#rownames(genexp1)<-genexp1$Description
rownames(genexp1)<-genexp1$Name


#Cambiamos los nombres de columumna por los nombres de genes
snames<-colnames(genexp1[,-1])
genexp2<-genexp1[,snames]



#####################
### Normalizacion ###
#####################

#Generamos la matriz de diseño
y<-DGEList(genexp2)
str(y)

#Normalizamos
y_norm<-calcNormFactors(y, method = "TMM")
y_norm1 <- cpm(y_norm, log=T)
#dim(y_norm1)

y_norm1<- y_norm1[1:100,1:200]

#################################
### Actividad de señalización ###
#################################

# https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/hipathia/inst/doc/hipathia-vignette.pdf

#Los nombres de las filas deben ser los ID de Entrez. Función translate_data para transformar IDs de genes en IDs de Entrez
trans_data <- translate_data(y_norm1, "hsa")

#la matriz de datos de expresión debe escalarse entre 0 y 1 antes de calcular los valores de activación de las sub-rutas
exp_data <- normalize_data(trans_data)

#boxplot(trans_data)
#boxplot(exp_data)


#carga 146 rutas
pathways <- load_pathways(species = "hsa")

#calcula el nivel de activación de las sub-rutas
results <- hipathia(exp_data, pathways, decompose = FALSE, verbose=FALSE)

#extraer el objeto de los valores de actividad de la señal 
path_vals <- get_paths_data(results)

#normalizar para comparar valores de señal entre sub-rutas
norm_paths<-normalize_paths(path_vals, pathways)
str(norm_paths)


















