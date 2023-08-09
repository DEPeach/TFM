##########################################################################################

# Python Script TFM - UEM
# Group_JMI_Rand

#Alumno: Sandra S. L. Roldán Pinzón
#Tutor: Carlos Loucera

# Implementación de un método basado en la entropía para la selección de variables 
#en modelos de aprendizaje máquina multi-respuesta para el reposicionamiento de fármacos.

#Input data: 
#           - Matriz de expresion genica de GTEx V8 normalizada y filtrada 
#           - KDTs del DrugBank Version 5.1.10

# Inputs (ORIGINAL) ## TO BE EDITED!!!!!!!!
#    X_data: n x d matrix X, with categorical values for n examples and d features
#    Y_labels: n x q matrix with the labels
#    topK: Number of features to be selected
#    distance: the distance measure that will be used for clustering the output space,
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Descripcion: 
#           - Implementacion del algoritmo de selección de características multivariante 
#             basado en entropía usando como guía la API de scikit-learn.
#           - Programa Group_JMI_Rand: Group_JMI_Rand algorithm for feature selection in multi-target problems

#Fuente: 
#        - https://www.mdpi.com/1099-4300/21/9/855
#        - https://github.com/sechidis/2019-Entropy-Multi-target-feature-selection/blob/master/Group_JMI_Rand.m

##########################################################################################


import numpy as np
import pandas as pd
import scipy as sp
import sklearn 



def Group_JMI_Rand(X_data,Y_labels, topK, distance):

    num_features = X_data.shape[1]
    num_labels = Y_labels.shape[1]
    num_ensemble = num_labels

    # Create the cluster ensebmle, as a first approach I will generate the same
    # number of targets as the initial, numLabels
    for index_label in range(num_ensemble):
        Y_labels_new(:,index_label) = kmedoids(Y_labels(:, datasample(1:num_labels,randi([ceil(num_labels/4) floor(3*num_labels/4)]),'Replace',false)),randi([4 16]), 'Distance', distance);
    


    score_per_feature = np.zeros(1,num_features)
    for index_feature in range(num_features):
        for index_label in range(num_ensemble):
            score_per_feature(index_feature) = score_per_feature(index_feature) + mi(X_data(:,index_feature),Y_labels_new(:,index_label))/sqrt(h(X_data(:,index_feature)) *h(Y_labels_new(:,index_label)));

    [val_max,selectedFeatures(1)]= max(score_per_feature)
    not_selected_features = setdiff(1:num_features,selectedFeatures)

    #cEfficient implementation of the second step, at this point I will store
    # the score of each feature. Whenever I select a feature I put NaN score
    score_per_feature = zeros(1,num_features)
    score_per_feature(selectedFeatures(1)) = NaN
    count = 2;

    while count<=topK:
        for index_feature_ns in len(not_selected_features):
            for index_label in range(num_ensemble): 
                score_per_feature(not_selected_features(index_feature_ns)) = score_per_feature(not_selected_features(index_feature_ns))+mi([X_data(:,not_selected_features(index_feature_ns)),X_data(:, selectedFeatures(count-1))], Y_labels_new(:,index_label))/sqrt(h([X_data(:,not_selected_features(index_feature_ns)),X_data(:, selectedFeatures(count-1))])*h(Y_labels_new(:,index_label)));
           
        [val_max,selectedFeatures(count)]= nanmax(score_per_feature)
        score_per_feature(selectedFeatures(count)) = NaN
        not_selected_features = setdiff(1:num_features,selectedFeatures)
        count = count+1
