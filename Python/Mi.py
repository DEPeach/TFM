##########################################################################################

# Python Script TFM - UEM
# Mi

#Alumno: Sandra S. L. Roldán Pinzón
#Tutor: Carlos Loucera

# Implementación de un método basado en la entropía para la selección de variables 
#en modelos de aprendizaje máquina multi-respuesta para el reposicionamiento de fármacos.

#Input data: 
#           - Matriz de expresion genica de GTEx V8 normalizada y filtrada 
#           - KDTs del DrugBank Version 5.1.10

#Descripcion: 
#           - Implementacion del algoritmo de selección de características multivariante 
#             basado en entropía usando como guía la API de scikit-learn.
#           - Programa MI: Estimate mutual information I(X;Y) between two categorical variables X,Y
#             X,Y can be matrices which are converted into a joint variable before computation

#Fuente: 
#        - https://www.mdpi.com/1099-4300/21/9/855
#        - https://github.com/sechidis/2019-Entropy-Multi-target-feature-selection/blob/master/mi.m

##########################################################################################

import numpy as np
import pandas as pd
import scipy as sp


X = [[2, 4, 4, 1],
      [6, 0, 0, 0],
      [6, 0, 0, 0]]

Y = [[1, 4, 4, 0],
      [5, 5, 5, 1],
      [5, 5, 5, 1]]

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
    p_X_p_Y=np.dot(p_XY.sum(axis=1), p_XY.sum(axis=0))


    #A=np.argwhere(p_XY)
    #B=np.argwhere(p_X_p_Y)
    #id_non_zero=np.where((A[:, None] == B).all(-1).any(1) == True)
    #MI = sum(sum( p_XY(id_non_zero) .* log(p_XY(id_non_zero) ./  p_X_p_Y(id_non_zero))))
    
    ##return MI