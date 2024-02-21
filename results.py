import glob
import re
import os
import pandas as pd

prefix = "PAX"

def get_model_data(filename):
    with open(filename) as f:
        lines = f.readlines()
        name = f.name.split("model_results.txt")[0]
        num_genes = 0
        if os.path.exists(name + "rna_paired_count_dummy_70.csv"):
            df = pd.read_csv(name + "rna_paired_count_dummy_70.csv")
            num_genes = df.shape[0]
        return {
            "name": name.split(f"/{prefix}_")[-1][0:-1],
            "r2": float(lines[3].split(":  ")[1]),
            "pearson": float(lines[4].split(":  ")[1]),
            "spearman": float(lines[5].split(":  ")[1]),
            "time": float(lines[6].split(":  ")[1]),
            "num_genes": num_genes
        }

df = pd.DataFrame(map(get_model_data, glob.glob(f"./trials/{prefix}_*/model_results.txt")))
df.sort_values("r2", ascending=False).to_csv(f"{prefix}_experiments.csv")
