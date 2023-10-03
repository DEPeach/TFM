import pandas as pd
import numpy as np
from disc_dataset_equalwidth import disc_dataset_equalwidth
from Group_JMI_Rand import Group_JMI_Rand
from Estabilidad import nogueria_test


np.random.seed(1026) 

#import data
X_org=pd.read_feather(path="/Users/macondo/Documents/GitHub/TFM/Raw_Data/pathvals_gtex-v8_edger-v3-40-0_hipathia-norm-v2-14-0.feather")
labels_org=pd.read_feather(path="/Users/macondo/Documents/GitHub/TFM/Raw_Data/labels.feater")

X=X_org.to_numpy()
Y_labels=labels_org.to_numpy()

X=X[:,1:733]
Y_labels=Y_labels[0:,1:26]

bins=5
X_inputs_disc = disc_dataset_equalwidth(X, bins)

topK = 10
Selected_with_Group_JMI_1 = Group_JMI_Rand(X_inputs_disc,Y_labels, topK, 'euclidean')
Selected_with_Group_JMI_2 = Group_JMI_Rand(X_inputs_disc,Y_labels, topK, 'hamming')



def estabilidad(n_samples, X, Y_labels, bins, topK):
    estab_matriz = np.zeros((n_samples, X.shape[1]))
    for sample in range(n_samples):
        random1=np.random.choice(X.shape[0], 50, replace=False)
        X_rd=X[random1,:]
        Y_labels_rd=Y_labels[random1,:]
        X_inputs_disc = disc_dataset_equalwidth(X_rd, bins)
        Selected_with_Group_JMI_1 = Group_JMI_Rand(X_inputs_disc,Y_labels_rd, topK, 'euclidean')
        estab_matriz[sample, Selected_with_Group_JMI_1] = 1
    return estab_matriz

estab_matriz=estabilidad(5, X, Y_labels, 5, 10)

nogueria_score=nogueria_test(estab_matriz)