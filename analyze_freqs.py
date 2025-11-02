import pandas as pd
import numpy as np

aa = pd.read_csv("results/aa_freqs.tsv", sep="\t", index_col=0)
di = pd.read_csv("results/dipep_freqs.tsv", sep="\t", index_col=0)

def top_variants(df, topn=10):
    var = df.var(axis=0)
    top = var.sort_values(ascending=False).head(topn)
    return top

print("Top 10 most variable amino acids:")
print(top_variants(aa))

print("\n Top 10 most variable dipeptides:")
print(top_variants(di))

