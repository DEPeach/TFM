##########################################################################################

# R Script TFM - UEM

#Alumno: Sandra S. L. Roldán Pinzón
#Tutor: Carlos Loucera

# Implementación de un método basado en la entropía para la selección de variables 
#en modelos de aprendizaje máquina multi-respuesta para el reposicionamiento de fármacos.

#Raw data: 
#         - GTEx V8
#         - DrugBank Version 5.1.10

#Descripcion: 
#           - Normalizacion de los perfiles de GTEx utilizando egdeR
#           - Calculo de los perfiles de actividad de señalización utilizando Hipathia
#           - Filtracion de los targets de drogas conocidos de la base de datos de DRUGBANK 
#             (Known Drug Targets, KDTs)

##########################################################################################

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
library(dplyr)


################################################################
####  Importacion de archivos y filtracion/cruzado de datos ####
################################################################

#Matriz de datos de expresion:
genexp <- read.delim(file="./Raw_Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2)

#Dimensiones de las matrices
dim(genexp)
genexp1<-genexp[,-2]   #elimina la columna description

#Comprobar si hay duplicados
sum(duplicated(genexp1$Name))


#Cambiamos los nombres de fila por los nombres de muestra
#variable name: nomenclatura de genes "Ensembl"
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
pathways <- hipathia::load_pathways("hsa")

#calcula el nivel de activación de las sub-rutas
results <- hipathia(exp_data, pathways, decompose = FALSE, verbose=FALSE)

#extraer el objeto de los valores de actividad de la señal 
path_vals <- get_paths_data(results, matrix = TRUE)

#normalizar para comparar valores de señal entre sub-rutas
norm_paths<-normalize_paths(path_vals, pathways)
str(norm_paths)




#################################################################################################################
#### Filtracion de los targets de drogas conocidos de la base de datos de DRUGBANK (Known Drug Targets, KDTs) ###
#################################################################################################################

# Carga Drugbank data
install.packages("dbparser")
library("dbparser")

drugbank_import <- parseDrugBank(db_path            = "./Raw_Data/drugbank_all_full_database_v5.1.8.xml",
                           drug_options       = c("groups"),
                           parse_salts        = FALSE,
                           parse_products     = FALSE,
                           references_options = NULL,
                           cett_options       = c("targets"))

## Carga drug groups data
drug_group <- drugbank_import$drugs$groups
drug_group$drugbank_id<-drug_group$`drugbank-id`

## Carga drug targets actions data
drug_targets <- drugbank_import$cett$targets$general_information
drug_targets <- drug_targets %>% rename(drugbank_id = parent_key)

## Carga drug targets uniprot data
drug_uniprot <- drugbank_import$cett$targets$polypeptides$external_identy
drug_uniprot <- drug_uniprot[drug_uniprot$resource=="UniProtKB",]
drug_uniprot <- drug_uniprot %>% rename(id = parent_key)
drug_uniprot <- drug_uniprot %>% rename(uniprot_ids = identifier)


## merge de los targest datasets
drug_targets2 <- merge(drug_targets,drug_uniprot)

## merge de los drug groups con los targes para crear la base de datos con la cual se van a filtrar los KDTs
drugbank_all <- merge(drug_group,drug_targets2, by="drugbank_id")


##  Total KDTs
distinct(drugbank_all[, c(1,10)]) %>% dim(.) ## 18467
length(unique(drugbank_all$drugbank_id)) ## 7540 drugs
length(unique(drugbank_all$uniprot_id)) ## 4717 drugs


# filtro: group= aprovado ó withdrawn + accion conocida = si + organismo = humanos. Total 2992 records
drugbank_all2 <- drugbank_all[-(base::grep("withdrawn",drugbank_all$group)),] %>% .[(base::grep("^approved|approved,", .$group)),] %>% 
  .[(base::grep("yes", .$known_action)),] %>% .[.$organism == "Humans",] 

dim(drugbank_all2) ## 2992  KDT-drug 
length(unique(drugbank_all2$uniprot_id)) ## 744 Uniprot IDs  




## Load genes stable translator used to Unipro IDs
entrez_uniprot <- read.delim(file = file.path(data_folder2, "genes_drugbank-v050108_mygene-20230120.tsv")) %>% .[.$uniprot_id %in% drugbank_app_action$uniprot_id , ]
dim(entrez_uniprot)

drugbank_app_action <- merge(drugbank_app_action, entrez_uniprot,) %>% .[-which(is.na(.$entrez_id)),] 
length(unique(drugbank_app_action$entrez_id)) 

saveRDS(drugbank_app_action, file.path(rds_folder, "drugbank_app_action_entrez.rds"))

drugbank_app_action_genes <- merge(drugbank_app_action, genes_tr, by.x = "entrez_id", by.y = "entrez" ) ## Add the symbol column from the stable translate table
length(unique(drugbank_app_action_genes$entrez_id)) ## We have all included KDTs 711

write.xlsx(drugbank_app_action_genes, file =  file.path(tables_folder, "supp_tabl3_drugbank518_filtered.xlsx"))

alldrug_byaction <- drugbank_app_action_genes[,c("name", "actions", "entrez_id")]
colnames(alldrug_byaction) <- c("drug", "drug_action", "KDT")  

















