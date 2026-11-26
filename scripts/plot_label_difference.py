"""Figure B: benign vs pathogenic means and GC-adjusted differences with CIs.

Run: python plot_label_difference.py [summary.csv] [output prefix]
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, "/mnt/disk2/vishnu/dna_tokenizers/scripts/tokenizer_analysis")
from common import LEARNED, COLORS

infile = sys.argv[1] if len(sys.argv) > 1 else "clinvar_stability_summary.csv"
prefix = sys.argv[2] if len(sys.argv) > 2 else "clinvar_fig_b"
LABELS = {"BPE": "Standard BPE"}
res = pd.read_csv(infile)
names = [LABELS.get(n, n) for n in LEARNED]
GROUPS = [("missense", "Missense variants"), ("intron", "Intronic variants")]

# Figure 1: averages, benign next to pathogenic
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
w = 0.38
for ax, (g, gtitle) in zip(axes, GROUPS):
    s = res[(res.group == g) & (res.metric == "bases_retok")].set_index("tokenizer")
    for i, name in enumerate(LEARNED):
        ax.barh(i - w / 2, s.loc[name, "mean_benign"], height=w, color="#7f8c9a")
        ax.barh(i + w / 2, s.loc[name, "mean_pathogenic"], height=w, color="#d1495b")
    ax.set_yticks(range(len(LEARNED)))
    ax.set_yticklabels(names)
    ax.set_ylim(len(LEARNED) - 0.5, -0.5)
    ax.set_title(gtitle, loc="left", fontweight="bold")
    ax.set_xlabel("bases re-tokenized per substitution (mean)")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].legend(handles=[Patch(color="#7f8c9a", label="benign"), Patch(color="#d1495b", label="pathogenic")], frameon=False)
fig.suptitle("Letters re-tokenized: benign vs pathogenic ClinVar variants", fontweight="bold")
fig.text(0.5, 0.01, "Means only; confidence ranges are on the difference figure. Fixed-length tokenizers are identical for both labels and not shown.",
         ha="center", fontsize=8, color="gray")
fig.tight_layout(rect=[0, 0.04, 1, 0.94])
fig.savefig(prefix + "_means.png", dpi=200)
fig.savefig(prefix + "_means.pdf")

# Figure 2: GC-adjusted difference with confidence range
metrics = [("bases_retok", "letters re-tokenized"), ("count_changed", "token count changes (fraction)")]
fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharey=True)
for r, (metric, mlabel) in enumerate(metrics):
    for c, (g, gtitle) in enumerate(GROUPS):
        ax = axes[r, c]
        s = res[(res.group == g) & (res.metric == metric)].set_index("tokenizer")
        for i, name in enumerate(LEARNED):
            m, lo, hi = s.loc[name, "gc_adjusted_diff"], s.loc[name, "ci_low"], s.loc[name, "ci_high"]
            crosses = lo <= 0 <= hi
            ax.barh(i, m, color=COLORS[name], alpha=0.4 if crosses else 1.0,
                    xerr=[[m - lo], [hi - m]], capsize=3, error_kw={"lw": 1})
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(range(len(LEARNED)))
        ax.set_yticklabels(names)
        ax.set_ylim(len(LEARNED) - 0.5, -0.5)
        ax.set_title(f"{gtitle}: {mlabel}", loc="left", fontweight="bold", fontsize=10)
        ax.set_xlabel("pathogenic minus benign (GC-adjusted)")
        ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Do tokenizers re-cut pathogenic variants more or less than benign ones?", fontweight="bold")
fig.text(0.5, 0.01, "Negative = pathogenic variants change less. Error bars: 95% range, resampled by gene. "
         "Faded bars: range includes zero. Fixed-length tokenizers are exactly zero and not shown.",
         ha="center", fontsize=8, color="gray")
fig.tight_layout(rect=[0, 0.04, 1, 0.95])
fig.savefig(prefix + "_difference.png", dpi=200)
fig.savefig(prefix + "_difference.pdf")
print("saved figures with prefix", prefix)
