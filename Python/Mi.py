##########################################################################################

# Python Script TFM - UEM
# Mi

#Alumno: Sandra S. L. Roldán Pinzón
#Tutor: Carlos Loucera

# Implementación de un método basado en la entropía para la selección de variables 
#en modelos de aprendizaje máquina multi-respuesta para el reposicionamiento de fármacos.

#Descripcion: 
#           - Implementacion del algoritmo de selección de características multivariante 
#             basado en entropía usando como guía la API de scikit-learn.
#           - Programa MI: Estima la información mutua I(X;Y) entre dos variables categóricas X,Y

#Fuente: 
#        - https://www.mdpi.com/1099-4300/21/9/855
#        - https://github.com/sechidis/2019-Entropy-Multi-target-feature-selection/blob/master/mi.m
##########################################################################################

import numpy as np
import pandas as pd
import scipy as sp


def mi(X,Y):
    C, ia, X1 = np.unique(X, return_index=True, return_inverse=True, axis=0)
    C, ia, Y1 = np.unique(Y, return_index=True, return_inverse=True, axis=0)
     
    C, ia, X2 = np.unique(X1, return_index=True, return_inverse=True)
    arity_X=len(ia)
    C, ia, Y2 = np.unique(Y1, return_index=True, return_inverse=True)
    arity_Y=len(ia)
    n = len(Y2)
    
    H, egdes = np.histogramdd(np.stack((X2, Y2), axis=1), bins = arity_Y)
    p_XY = H/n
    p_X_p_Y=np.transpose(np.array([p_XY.sum(axis=1)]))@np.array([p_XY.sum(axis=0)])
    arg_XY=np.argwhere(p_XY)
    arg_X_Y=np.argwhere(p_X_p_Y)
    mask=(arg_XY[:, None] == arg_X_Y).all(-1).any(1)
    id_non_zero=arg_XY[mask]
    
    MI = sum(np.multiply(p_XY[tuple(zip(*id_non_zero))],np.log(p_XY[tuple(zip(*id_non_zero))] / p_X_p_Y[tuple(zip(*id_non_zero))])))
        
    return MI