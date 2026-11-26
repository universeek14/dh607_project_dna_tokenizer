"""Keep at most 10 variants per (group, label, gene), seed 42 (10,714 of 12,000 kept).

Input: clinvar_stability_long.tsv. Output: clinvar_stability_long_capped.tsv.
"""
import pandas as pd

d = pd.read_csv("clinvar_stability_long.tsv", sep="\t")
d["gene"] = d["gene"].fillna("NA")
base = d[d.tokenizer == "char"][["idx", "group", "label", "gene"]]
base = base.sample(frac=1, random_state=42)
keep = base.groupby(["group", "label", "gene"]).head(10)["idx"]
d[d.idx.isin(keep)].to_csv("clinvar_stability_long_capped.tsv", sep="\t", index=False)
print(len(keep), "variants kept")
print(base[base.idx.isin(keep)].groupby(["group", "label"]).size())
