"""Pathogenic - benign stability differences for the 9 pretraining tokenizers.

Same method as summarize_stability.py. Output: clinvar_stability_summary_ablation.csv.
"""
import numpy as np
import pandas as pd

d = pd.read_csv("clinvar_stability_long_ablation.tsv", sep="\t")
d["count_changed"] = (d["n_ref"] != d["n_alt"]).astype(int)
d["gene"] = d["gene"].fillna("NA")
rng = np.random.default_rng(42)
B = 500


def adj_diff(y, lab, gc):
    X = np.column_stack([np.ones(len(y)), lab, gc - gc.mean()])
    return np.linalg.lstsq(X, y, rcond=None)[0][1]


out = []
for name in d["tokenizer"].unique():
    for g in ["missense", "intron"]:
        s = d[(d.tokenizer == name) & (d.group == g)].reset_index(drop=True)
        lab = s["label"].values.astype(float)
        gc = s["gc"].values
        genes = s["gene"].values
        uniq = np.unique(genes)
        idx = [np.where(genes == u)[0] for u in uniq]
        for metric in ["bases_retok", "count_changed"]:
            y = s[metric].values.astype(float)
            boots = []
            for _ in range(B):
                pick = rng.integers(0, len(uniq), len(uniq))
                ii = np.concatenate([idx[k] for k in pick])
                boots.append(adj_diff(y[ii], lab[ii], gc[ii]))
            out.append({
                "tokenizer": name, "group": g, "metric": metric,
                "mean_benign": y[lab == 0].mean(), "mean_pathogenic": y[lab == 1].mean(),
                "raw_diff": y[lab == 1].mean() - y[lab == 0].mean(),
                "gc_adjusted_diff": adj_diff(y, lab, gc),
                "ci_low": np.percentile(boots, 2.5), "ci_high": np.percentile(boots, 97.5),
            })

res = pd.DataFrame(out)
res.to_csv("clinvar_stability_summary_ablation.csv", index=False)
print(res.round(3).to_string(index=False))
