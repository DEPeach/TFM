##########################################################################################

# Python Script TFM - UEM
# Group_JMI_Rand

#Alumno: Sandra S. L. Roldán Pinzón
#Tutor: Carlos Loucera

# Implementación de un método basado en la entropía para la selección de variables 
#en modelos de aprendizaje máquina multi-respuesta para el reposicionamiento de fármacos.

# Inputs:
#           - Datos_X: matriz con valores categóricos para n ejemplos y d características.
#           - Y_labels: matriz con las etiquetas
#           - topK: número de características a seleccionar
#           - distancia: la medida de distancia que se utilizará para agrupar el espacio de salida

#Descripcion: 
#           - Implementacion del algoritmo de selección de características multivariante 
#             basado en entropía usando como guía la API de scikit-learn.
#           - Programa Group_JMI_Rand: Algoritmo Group_JMI_Rand para la selección de características en problemas multiobjetivo

#Fuente: 
#        - https://www.mdpi.com/1099-4300/21/9/855
#        - https://github.com/sechidis/2019-Entropy-Multi-target-feature-selection/blob/master/Group_JMI_Rand.m
##########################################################################################


import numpy as np
import pandas as pd
import scipy as sp
import math
import random
import sklearn 
from sklearn_extra.cluster import KMedoids
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics import pairwise_distances
from scipy.stats import entropy
import Mi

def Group_JMI_Rand(X_data,Y_labels, topK, distance):
    
    num_features = X_data.shape[1]
    num_labels = Y_labels.shape[1]
    num_ensemble = num_labels
    Y_labels_new = np.zeros_like(Y_labels)
    
    for index_label in range(num_ensemble):
        #random1 = np.random.randint(math.ceil(num_labels/4),math.floor(3*num_labels/4)+1)
        random1 = np.random.choice(np.arange(num_labels), size=np.random.randint(np.ceil(num_labels/4), np.floor(3*num_labels/4)), replace=False)
        random2 = np.random.randint(4,16+1)
        #random_sample=random.sample([x for x in range(0, num_labels, 1)],random1)
        random_sample = Y_labels[:, random1]
        #Y_labels_new = KMedoids(n_clusters=random2, random_state=0, method="pam").fit(np.array([row[random_sample] for row in Y_labels]))
        kmedoids = KMedoids(n_clusters=random2, random_state=None, metric=distance).fit(random_sample)
        Y_labels_new[:, index_label] = kmedoids.labels_
        
    score_per_feature = np.zeros((1,num_features),dtype=float)
    
    for index_feature in range(num_features):
        for index_label in range(num_ensemble):
            # Entropia de X_data y Y_labels_new
            h_X = entropy(X_data[:, index_feature])
            h_Y = entropy(Y_labels_new[:, index_label])
            score_per_feature[0, index_feature] += Mi.mi(X_data[:, index_feature],Y_labels_new[:, index_label])/np.sqrt(h_X * h_Y)
                        
    valid_indices = np.where(~np.isnan(score_per_feature))[1]
    selectedFeatureIndex = valid_indices[np.argmax(score_per_feature[0,valid_indices])]
    val_max = score_per_feature[0,selectedFeatureIndex]
    all_features = set(range(0, num_features))
    selectedFeatures = [selectedFeatureIndex] 
    not_selected_features = sorted(list(all_features - set(selectedFeatures)))
    score_per_feature = np.zeros((1,num_features))
    score_per_feature[0,selectedFeatures[0]] = np.nan
    count = 2
    
    while count<=topK:
        for index_feature_ns in range(len(not_selected_features)):
            for index_label in range(num_ensemble):       
                combined_features = np.stack((X_data[:, not_selected_features[index_feature_ns]], X_data[:, selectedFeatures[count - 2]]), axis=1)
                combined_features_ent = np.hstack((X_data[:, not_selected_features[index_feature_ns]], X_data[:, selectedFeatures[count - 2]]))
                h_X = entropy(combined_features_ent)
                h_Y = entropy(Y_labels_new[:, index_label])
                score_per_feature[0, not_selected_features[index_feature_ns]] += Mi.mi(combined_features,Y_labels_new[:, index_label]) / np.sqrt(h_X * h_Y)                       
        valid_indices = np.where(~np.isnan(score_per_feature))[1]
        selectedFeatureIndex = valid_indices[np.argmax(score_per_feature[0,valid_indices])]
        val_max = score_per_feature[0,selectedFeatureIndex]
        score_per_feature[0,selectedFeatureIndex] = np.nan
        selectedFeatures.append(selectedFeatureIndex)  
        not_selected_features = sorted(list(set(range(0, num_features)) - set(selectedFeatures)))
        count += 1
        
    return selectedFeatures


