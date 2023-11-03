import pathlib
import joblib
import numpy as np
import pandas as pd
from drexml.datasets import load_physiological_circuits
from joblib import Parallel, delayed
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split

from disc_dataset_equalwidth import disc_dataset_equalwidth
from Group_JMI_Rand import Group_JMI_Rand
def get_stab_results(n_samples, features, targets, bins, topK):
    np.random.seed(0)
    def run_one(sample):
        X_train, X_test, Y_train, Y_test = train_test_split(
            features, targets, test_size=0.5, random_state=sample
        )
        np.random.seed(sample)
        X_train_disc = disc_dataset_equalwidth(X_train.to_numpy(), bins)
        Selected_with_Group_JMI_1 = Group_JMI_Rand(
            X_train_disc, Y_train.to_numpy(), topK, "euclidean"
        )
        estab_matriz = pd.Series(np.zeros(X_train.shape[1]),
index=X_train.columns)
        estab_matriz.iloc[Selected_with_Group_JMI_1] = 1
        estab_matriz.name = f"seed_{sample}"
        X_train_filtered = X_train.iloc[:, Selected_with_Group_JMI_1]
        X_test_filtered = X_test.iloc[:, Selected_with_Group_JMI_1]
        model = RidgeCV()
        model.fit(X_train_filtered, Y_train)
        Y_test_hat = model.predict(X_test_filtered)
        coef = pd.DataFrame(
            model.coef_.astype(bool),
            index=Y_train.columns,
            columns=X_train_filtered.columns,
)
        r2_vals = r2_score(Y_test, Y_test_hat, multioutput="raw_values")
        r2_vals = pd.Series(r2_vals, index=Y_test.columns)
        r2_vals.name = f"seed_{sample}"
        return estab_matriz, r2_vals, coef
    res = Parallel(n_jobs=-1)(delayed(run_one)(i) for i in range(n_samples))
    return res

if __name__ == "__main__":
    project_path = pathlib.Path(".", "..", "Raw_Data")
    data_path = pathlib.Path("~/.data/zenodo/6020480/20230904-fix01")
    X = pd.read_feather(
        path=data_path.joinpath("gexp_gtex-v8_edger-v3-40-0.feather")
    ).set_index("index")
    X.columns = X.columns.str.replace("X", "")
    kdt_entrez_list = (
        pd.read_csv(
            project_path.joinpath("Drugbank_KDTs_mygene.tsv"), sep="\t",
dtype=str )
        .entrez_id.dropna()
.unique() )
    entrez_to_symbol = (
        pd.read_csv(
            project_path.joinpath("gexp_gtex-v8_edger-v3-40-0_mygene.tsv"),
            sep="\t",
            dtype="str",
        )
        .set_index("entrez_id")
        .loc[:, "symbol_id"]
        .to_dict()
)
    X = X[X.columns.intersection(kdt_entrez_list)]
    X = X.rename(columns=entrez_to_symbol)
    physio_circuits = load_physiological_circuits()
    disease_map =
pd.read_feather(project_path.joinpath("circuits_AF.feater"))
    disease_map["Hipathia_code"] = (
        disease_map["Hipathia_code"].str.replace("-", ".").str.replace(" ",
".") )
disease_map = disease_map.set_index("Hipathia_code").loc[:, "Circuit"] # disease_map = disease_map.loc[disease_map.index.intersection(physio_circuits)]
Y = pd.read_feather(
    path=data_path.joinpath(
        "pathvals_gtex-v8_edger-v3-40-0_hipathia-v2-14-0.feather"
    )
).set_index("index")
Y = Y[Y.columns.intersection(disease_map.index)]
Y = Y.rename(columns=disease_map)
res = get_stab_results(100, features=X, targets=Y, bins=10, topK=k_top)
joblib.dump(res, "resultados.pkl")