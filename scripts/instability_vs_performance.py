"""Analysis 4 (exploratory): SNV instability vs downstream score, 9 tokenizers.

Performance values (seed 1) are typed in from pretrain_small_bert/README.md.
Spearman rho over 9 and over the 7 BPE-family tokenizers. n is small; see README.md.
Output: instability_vs_performance.csv/png.
"""
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Seed-1 numbers copied from the README table of pretrain_small_bert
perf = pd.DataFrame({
    "tokenizer":     ["char", "kmer6", "bpe", "len2_only", "strat_only", "evolen",
                      "evolen_rebuilt", "random_strat", "func_strat"],
    "clinvar_probe": [0.661, 0.662, 0.630, 0.628, 0.631, 0.617, 0.609, 0.629, 0.622],
    "gue_splice":    [82.7, 80.8, 73.3, 68.8, 72.5, 71.2, 67.9, 69.3, 68.9],
    "gue_mean":      [79.7, 76.6, 72.9, 70.9, 73.2, 72.5, 70.6, 72.0, 71.8],
}).set_index("tokenizer")

d = pd.read_csv("clinvar_stability_long_ablation.tsv", sep="\t")
d["changed"] = (d.n_ref != d.n_alt).astype(float)
inst = d.groupby("tokenizer").agg(retok=("bases_retok", "mean"), count_changed=("changed", "mean"))
t = perf.join(inst)
t.to_csv("instability_vs_performance.csv")
print(t.round(3))

family = [x for x in t.index if x not in ("char", "kmer6")]
for label, sub in [("all nine", t), ("BPE family only (7)", t.loc[family])]:
    print("\n" + label)
    for m in ["clinvar_probe", "gue_splice", "gue_mean"]:
        for x in ["retok", "count_changed"]:
            rho, p = spearmanr(sub[x], sub[m])
            print(f"  {x:14s} vs {m:14s} rho {rho:+.2f} (p {p:.2f})")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, (m, lab) in zip(axes, [("gue_splice", "GUE splice MCC"), ("clinvar_probe", "ClinVar probe AUROC")]):
    ax.scatter(t["retok"], t[m], color="#1e88e5")
    for name, r in t.iterrows():
        ax.annotate(name, (r["retok"], r[m]), fontsize=8, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("letters re-tokenized per ClinVar substitution (mean)")
    ax.set_ylabel(lab)
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Tokenizer instability vs downstream score (9 tokenizers, seed 1)", fontweight="bold")
fig.text(0.5, 0.005, "Exploratory: nine points, resolution and vocabulary differ together.", ha="center", fontsize=8, color="gray")
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.savefig("instability_vs_performance.png", dpi=200)
