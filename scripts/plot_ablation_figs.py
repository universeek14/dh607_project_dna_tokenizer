import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORDER = ["char", "kmer6", "bpe", "len2_only", "strat_only", "evolen",
         "evolen_rebuilt", "random_strat", "func_strat"]
SHOW = {"char": "char (1 bp)", "kmer6": "6-mer", "bpe": "BPE",
        "len2_only": "BPE vocab + length decoding",
        "strat_only": "EvoLen vocab + frequency decoding",
        "evolen": "EvoLen", "evolen_rebuilt": "EvoLen rebuilt (phyloP strata)",
        "random_strat": "EvoLen, shuffled strata", "func_strat": "EvoLen, functional strata"}
COL = {"char": "#9e9e9e", "kmer6": "#9e9e9e", "bpe": "#1e88e5", "len2_only": "#fb8c00",
       "strat_only": "#8e24aa", "evolen": "#e53935", "evolen_rebuilt": "#7f0000",
       "random_strat": "#6d4c41", "func_strat": "#00897b"}

# Figure 1: letters re-tokenized, all variants, 95% range resampled by gene
d = pd.read_csv("clinvar_stability_long_ablation.tsv", sep="\t")
d["gene"] = d["gene"].fillna("NA")
rng = np.random.default_rng(42)
fig, ax = plt.subplots(figsize=(9, 5))
for i, name in enumerate(ORDER):
    g = d[d.tokenizer == name].groupby("gene")["bases_retok"].agg(["sum", "count"])
    s, c = g["sum"].values, g["count"].values
    boots = []
    for _ in range(1000):
        p = rng.integers(0, len(g), len(g))
        boots.append(s[p].sum() / c[p].sum())
    m = s.sum() / c.sum()
    lo, hi = np.percentile(boots, [2.5, 97.5])
    ax.barh(i, m, color=COL[name], xerr=[[m - lo], [hi - m]], capsize=2, error_kw={"lw": 0.8})
    ax.text(hi + 0.2, i, f"{m:.2f}", va="center", fontsize=9)
ax.set_yticks(range(len(ORDER)))
ax.set_yticklabels([SHOW[n] for n in ORDER])
ax.invert_yaxis()
ax.set_xlabel("bases re-tokenized per ClinVar substitution (mean)")
ax.set_title("Tokenizer response to real ClinVar substitutions: the nine pretraining tokenizers",
             loc="left", fontweight="bold", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
fig.text(0.5, 0.005, "Error bars: 95% range, resampled by gene. Measures tokenizer behaviour on the input, not model predictions.",
         ha="center", fontsize=8, color="gray")
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig("clinvar_fig_ablation_means.png", dpi=200)
fig.savefig("clinvar_fig_ablation_means.pdf")

# Figure 2: pathogenic minus benign, GC-adjusted, from the summary table
res = pd.read_csv("clinvar_stability_summary_ablation.csv")
names = [n for n in ORDER if n not in ("char", "kmer6")]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, (g, title) in zip(axes, [("missense", "Missense variants"), ("intron", "Intronic variants")]):
    s = res[(res.group == g) & (res.metric == "bases_retok")].set_index("tokenizer")
    for i, name in enumerate(names):
        m, lo, hi = s.loc[name, "gc_adjusted_diff"], s.loc[name, "ci_low"], s.loc[name, "ci_high"]
        ax.barh(i, m, color=COL[name], alpha=0.4 if lo <= 0 <= hi else 1.0,
                xerr=[[m - lo], [hi - m]], capsize=3, error_kw={"lw": 1})
        ax.text(hi + 0.05, i, f"{m:.2f}", va="center", fontsize=8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([SHOW[n] for n in names])
    ax.invert_yaxis()
    ax.set_title(title + ": letters re-tokenized", loc="left", fontweight="bold", fontsize=10)
    ax.set_xlabel("pathogenic minus benign (GC-adjusted)")
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Do the pretraining tokenizers re-cut pathogenic variants more or less than benign ones?", fontweight="bold")
fig.text(0.5, 0.005, "Negative = pathogenic variants change less. Error bars: 95% range, resampled by gene. Faded bars: range includes zero. Fixed-length tokenizers are exactly zero and not shown.",
         ha="center", fontsize=7, color="gray")
fig.tight_layout(rect=[0, 0.03, 1, 0.94])
fig.savefig("clinvar_fig_ablation_difference.png", dpi=200)
fig.savefig("clinvar_fig_ablation_difference.pdf")
print("saved")
