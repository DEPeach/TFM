##########################################################################################

# Python Script TFM - UEM
# Final_files

#Alumno: Sandra S. L. Roldán Pinzón
#Tutor: Carlos Loucera

# Implementación de un método basado en la entropía para la selección de variables 
#en modelos de aprendizaje máquina multi-respuesta para el reposicionamiento de fármacos.

#Input data: 
#           - Matriz de expresion genica de GTEx V8 normalizada y filtrada 
#           - KDTs del DrugBank Version 5.1.10

#Descripcion: 
#           - Preprocesamiento de los archivos finales para el analisis
#           - Conversión de nomenclatura de genes a través de my-gene -id

#Fuente: 
#        - https://github.com/babelomics/drexml-retinitis/blob/main/scripts/py/parser.py

##########################################################################################

import pandas as pd
from biothings_client import get_client


#import data



# definition of the function convert_gene_ids

def convert_gene_ids(gene_ids, source="entrezgene", target="uniprot,symbol"):
    """Myege wrapper."""
    renamer = {
        "uniprot": "uniprot_id",
        "entrezgene": "entrez_id",
        "symbol": "symbol_id",
    }
    client = get_client("gene")
    genes_converted = client.querymany(
        gene_ids,
        scopes=source,
        fields=target,
        species="human",
        as_dataframe=True,
    )

    genes_converted = genes_converted.reset_index(names=[source])
    cols_query = genes_converted.columns.isin(["uniprot", "entrezgene", "symbol"])
    genes_converted = genes_converted.loc[:, cols_query]
    genes_converted = genes_converted.rename(columns=renamer)

    return genes_converted



# definition of the function translate

def translate(path, output, kind):
    """Gene translation tool using mygene."""
    print("Running mygene translation tool.")

    path = Path(path)
    output = Path(output)

    kind = kind.lower()
    if kind == "drugbank":
        data = pd.read_csv(path, sep="\t")
        ids = data["uniprot_id"].unique()
        this_source = "uniprot"
        this_target = "entrezgene"
    elif kind == "gtex":
        data = pd.read_feather(path)
        ids = data.columns[data.columns.str.startswith("X")].str.replace("X", "")
        this_source = "entrezgene"
        this_target = "symbol"

    genes_df = convert_gene_ids(ids, source=this_source, target=this_target)
    genes_df.to_csv(output, sep="\t", index=False)
    print(f"Wrote {output}")
    