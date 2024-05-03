import gffutils
import pandas as pd

def annotate_pepper(db_fn):
    ann = pd.read_csv("2.0.renamed.annotation.final.tsv", sep="\t", index_col=[0])
    ann["Annotation"] = ann["Annotation"].astype(object)
    ann["Existence"] = ann["Existence"].astype(object)

    db = gffutils.FeatureDB(db_fn)
    cols = ann.columns

    new_entries = list()

    for gene in db.features_of_type("gene"):
        gene_ID = gene.attributes["ID"][0]
        _ = gene.attributes.pop("Name", None)

        if gene_ID in ann.index:
            data = ann.loc[gene_ID]

            for col in cols:
                val = data[col]
                #if not pd.isna(val):
                gene.attributes[col] = val if not pd.isna(val) else ""

        new_entries.append(gene)

    db = db.update(new_entries, merge_strategy="replace", make_backup=False)