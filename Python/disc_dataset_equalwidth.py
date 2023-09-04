##########################################################################################

# Python Script TFM - UEM
# Group_JMI_Rand

#Alumno: Sandra S. L. Roldán Pinzón
#Tutor: Carlos Loucera

# Implementación de un método basado en la entropía para la selección de variables 
#en modelos de aprendizaje máquina multi-respuesta para el reposicionamiento de fármacos.

# Inputs:
#           - X: matriz n x d X, con valores categóricos para n ejemplos y d características
#           - bins: el número de categorías

#Descripcion: 
#           - Implementacion del algoritmo de selección de características multivariante 
#             basado en entropía usando como guía la API de scikit-learn.
#           - Esta función discretiza cada característica de un conjunto de datos determinado de manera que tenga la misma anchura

#Fuente: 
#        - https://www.mdpi.com/1099-4300/21/9/855
#        - https://github.com/sechidis/2019-Entropy-Multi-target-feature-selection/blob/master/disc_dataset_equalwidth.m 

##########################################################################################


import numpy as np
import pandas as pd
import scipy as sp


def disc_dataset_equalwidth(X, bins):
    
    new_data = np.zeros_like(X)
    
    for fnum in range(X.shape[1]):
        
        if len(np.unique(X[:, fnum])) <= bins:
            _, _, new_data[:, fnum] = np.unique(X[:, fnum], return_inverse=True)
        else:
            feat = X[:, fnum]
            minval = np.min(feat)
            width = abs(np.max(feat) - minval) / bins

            boundaryend = np.zeros(bins)
            for n in range(bins):
                boundaryend[n] = minval +((n + 1) * width)
            boundaryend[-1] = boundaryend[-1] + 1
            
            lastboundaryend = minval
            newfeature = np.zeros_like(feat)
            
            for n in range(bins):
                indices = np.where((feat >= lastboundaryend) & (feat < boundaryend[n]))[0]
                newfeature[indices] = n + 1
                lastboundaryend = boundaryend[n]

            _, _, new_data[:, fnum] = np.unique(newfeature, return_inverse=True)
 