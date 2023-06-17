################################################

# Script TFM - UEM

#Alumno: Sandra Roldan 
###############################################

########################
### INSTALL PACKAGES ###
########################

if (!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager")}
if (!requireNamespace("edgeR", quietly=TRUE)){BiocManager::install("edgeR")}
if (!requireNamespace("tidyverse", quietly=TRUE)){install.packages("tidyverse")}


###########################
### PREPARE ENVIRONMENT ###
###########################

library(edgeR)
library(BiocManager)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(genefu)
library(tidyverse)


################################################################
####  Importacion de archivos y filtracion/cruzado de datos ####
################################################################

#Matriz de datos de expresion:
genexp <- read.delim(file="./Raw_Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2)



#Dimensiones de las matrices
dim(genexp)
genexp1<-genexp[,-1]   #elimina la columna name
genexp1<-genexp[,-2]   #elimina la columna description






#Comprobar si hay duplicados
sum(duplicated(genexp1$Description))
sum(duplicated(genexp1$Name))


#Cambiamos los nombres de fila por los nombres de muestra
rownames(genexp1)<-genexp1$Description
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
y_norm <- cpm(y_norm, log=T)
str(y_norm)
