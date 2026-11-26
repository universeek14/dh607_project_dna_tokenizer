"""Figure A: bases re-tokenized, token-count change, net split/merge per tokenizer.

Means with 95% gene-block bootstrap CIs (1,000 resamples).
Run: python plot_stability_deck_style.py [long.tsv] [output prefix]
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, "/mnt/disk2/vishnu/dna_tokenizers/scripts/tokenizer_analysis")
from common import TOKENIZERS, COLORS

infile = sys.argv[1] if len(sys.argv) > 1 else "clinvar_stability_long.tsv"
prefix = sys.argv[2] if len(sys.argv) > 2 else "clinvar_fig_a"
LABELS = {"char": "Character", "3mer": "3-mer", "6mer": "6-mer", "NT": "NT-v2", "BPE": "Standard BPE"}
FIXED = {"char", "3mer", "6mer", "NT"}
B = 1000

d = pd.read_csv(infile, sep="\t")
d["gene"] = d["gene"].fillna("NA")
d["changed"] = (d["n_ref"] != d["n_alt"]).astype(float)
d["split"] = (d["n_alt"] > d["n_ref"]).astype(float)
d["merge"] = (d["n_alt"] < d["n_ref"]).astype(float)
n_variants = d["idx"].nunique()
rng = np.random.default_rng(42)


def boot_ci(sums, counts):
    G = len(sums)
    vals = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, G, G)
        vals[b] = sums[pick].sum() / counts[pick].sum()
    return np.percentile(vals, [2.5, 97.5])


stats = {}
for name in TOKENIZERS:
    g = d[d.tokenizer == name].groupby("gene")
    cnt = g.size().values.astype(float)
    stats[name] = {}
    for m in ["bases_retok", "changed", "split", "merge"]:
        sm = g[m].sum().values
        lo, hi = boot_ci(sm, cnt)
        stats[name][m] = (sm.sum() / cnt.sum(), lo, hi)


def draw_simple(ax, metric, title, xlabel, show_names):
    xmax = max(stats[n][metric][2] for n in TOKENIZERS)
    for i, name in enumerate(TOKENIZERS):
        m, lo, hi = stats[name][metric]
        ax.barh(i, m, color=COLORS[name], xerr=[[m - lo], [hi - m]], capsize=2, error_kw={"lw": 0.8})
        ax.text(hi + 0.02 * xmax, i, f"{m:.2f}", va="center", fontsize=9)
    ax.set_yticks(range(len(TOKENIZERS)))
    ax.set_yticklabels([LABELS.get(n, n) for n in TOKENIZERS] if show_names else [])
    ax.invert_yaxis()
    ax.set_xlim(0, xmax * 1.2)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.spines[["top", "right"]].set_visible(False)


fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
draw_simple(axes[0], "bases_retok", "a  Bases re-tokenized", "bases re-tokenized per substitution (mean)", True)
draw_simple(axes[1], "changed", "b  Token-count change", "fraction of substitutions that change the number of tokens", False)

ax = axes[2]
xmax = max(max(stats[n]["split"][2], stats[n]["merge"][2]) for n in TOKENIZERS)
for i, name in enumerate(TOKENIZERS):
    if name in FIXED:
        ax.text(0.002, i, "0 / 0 (never changes count)", va="center", fontsize=8, color="gray")
        continue
    for off, key, col in [(-0.2, "split", "#4c72b0"), (0.2, "merge", "#dd8452")]:
        m, lo, hi = stats[name][key]
        ax.barh(i + off, m, height=0.38, color=col, xerr=[[m - lo], [hi - m]], capsize=2, error_kw={"lw": 0.8})
ax.set_yticks(range(len(TOKENIZERS)))
ax.set_yticklabels([])
ax.invert_yaxis()
ax.set_xlim(0, xmax * 1.2)
ax.set_title("c  Net split vs net merge", loc="left", fontweight="bold")
ax.set_xlabel("fraction of substitutions")
ax.legend(handles=[Patch(color="#4c72b0", label="net split (more tokens)"),
                   Patch(color="#dd8452", label="net merge (fewer tokens)")], frameon=False, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("Tokenizer response to real ClinVar substitutions", fontweight="bold")
fig.text(0.5, 0.01, f"{n_variants} ClinVar single-letter variants in 1-kb windows (missense and intronic, both labels). "
         "Error bars: 95% range, resampled by gene. Measures tokenizer behaviour on the input, not model predictions.",
         ha="center", fontsize=8, color="gray")
fig.tight_layout(rect=[0, 0.04, 1, 0.94])
fig.savefig(prefix + ".png", dpi=200)
fig.savefig(prefix + ".pdf")
print("saved", prefix + ".png and .pdf")
