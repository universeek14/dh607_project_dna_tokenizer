"""Analysis 2: analysis B (partial Spearman of tokens/kb) computed inside each region.

Features: mean phyloP, GC, repeat fraction, 3-mer entropy, each adjusted for the other
three. Data: data/windows.parquet (tokenizer-analysis window sample, 7 region classes).
CI: 10-Mb block bootstrap within region. Output: B_by_region/.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, "/mnt/disk2/vishnu/dna_tokenizers/scripts/tokenizer_analysis")
from common import (COLORS, LEARNED, REGION_ORDER, load_tokenizers, load_windows, log,
                    savefig, tokenize_windows)

ap = argparse.ArgumentParser()
ap.add_argument("--boot", type=int, default=300)
ap.add_argument("--out", default="B_by_region")
args = ap.parse_args()
OUT = Path(args.out)
OUT.mkdir(parents=True, exist_ok=True)

FEATURES = ["phylop_mean", "gc", "rep_frac", "kmer3_entropy"]
FLABEL = {"phylop_mean": "mean\nphyloP", "gc": "GC\nfraction",
          "rep_frac": "repeat\nfraction", "kmer3_entropy": "3-mer\nentropy"}
SHOW = {"BPE": "Standard BPE"}


def residual_rank(y, X):
    ry = rankdata(y)
    RX = np.column_stack([rankdata(c) for c in X.T] + [np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(RX, ry, rcond=None)
    return ry - RX @ beta


def partial_rho(y, X, j):
    others = [k for k in range(X.shape[1]) if k != j]
    a = residual_rank(y, X[:, others])
    b = residual_rank(X[:, j], X[:, others])
    return spearmanr(a, b)[0]


win = load_windows()
toks = load_tokenizers(LEARNED)
regions = [r for r in REGION_ORDER if r in set(win["region"])]
print("regions found:", regions)
rng = np.random.default_rng(42)
rows = []

for name, tok in toks.items():
    t = tokenize_windows(win, tok)
    w = win[["region", "chrom", "start"] + FEATURES].copy()
    w["tokens_per_kb"] = t.groupby("window").size().reindex(range(len(win))).values
    w = w.dropna().reset_index(drop=True)
    w["block"] = w["chrom"].astype(str) + "_" + (w["start"] // 10_000_000).astype(int).astype(str)
    subsets = [("All regions", w)] + [(r, w[w["region"] == r]) for r in regions]
    for label, s in subsets:
        s = s.reset_index(drop=True)
        y = s["tokens_per_kb"].values
        X = s[FEATURES].values
        groups = [np.asarray(ix) for ix in s.groupby("block").indices.values()]
        boots = np.empty((args.boot, len(FEATURES)))
        for b in range(args.boot):
            pick = rng.integers(0, len(groups), len(groups))
            ii = np.concatenate([groups[k] for k in pick])
            for j in range(len(FEATURES)):
                boots[b, j] = partial_rho(y[ii], X[ii], j)
        for j, f in enumerate(FEATURES):
            rows.append({"tokenizer": name, "region": label, "feature": f,
                         "partial_rho": partial_rho(y, X, j),
                         "ci_low": np.percentile(boots[:, j], 2.5),
                         "ci_high": np.percentile(boots[:, j], 97.5),
                         "n_windows": len(s), "n_blocks": len(groups)})
    log(f"{name} done")

res = pd.DataFrame(rows)
res.to_csv(OUT / "partial_correlations_by_region.csv", index=False)
print(res.round(3).to_string(index=False))

# Figure 1: panel a of the deck, one small panel per region
panels = ["All regions"] + regions
fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharey=True)
x = np.arange(len(FEATURES))
bw = 0.8 / len(toks)
for ax, label in zip(axes.flat, panels):
    for i, name in enumerate(toks):
        g = res[(res.tokenizer == name) & (res.region == label)].set_index("feature").loc[FEATURES]
        pos = x + (i - (len(toks) - 1) / 2) * bw
        ax.bar(pos, g["partial_rho"], bw, color=COLORS[name], label=SHOW.get(name, name),
               yerr=[g["partial_rho"] - g["ci_low"], g["ci_high"] - g["partial_rho"]],
               capsize=1.5, error_kw={"lw": 0.7})
    n = res[(res.region == label)]["n_windows"].iloc[0]
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([FLABEL[f] for f in FEATURES], fontsize=8)
    ax.set_title(f"{label} (n = {n})", loc="left", fontweight="bold", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
for ax in axes[:, 0]:
    ax.set_ylabel("partial Spearman rho with tokens/kb")
for ax in axes.flat[len(panels):]:
    ax.axis("off")
handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.97, 0.12), frameon=False)
fig.suptitle("What sets tokenization resolution? Each feature adjusted for the other three, by region",
             fontweight="bold")
fig.text(0.5, 0.01, "Error bars: 95% range, resampled over 10-Mb genome blocks within each region "
         f"({args.boot} resamples). Fixed tokenizers have constant resolution and are not shown.",
         ha="center", fontsize=8, color="gray")
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
savefig(fig, OUT / "fig_by_region_all_features")

# Figure 2: conservation only, regions side by side
fig, ax = plt.subplots(figsize=(11, 4.5))
x = np.arange(len(panels))
for i, name in enumerate(toks):
    g = res[(res.tokenizer == name) & (res.feature == "phylop_mean")].set_index("region").loc[panels]
    pos = x + (i - (len(toks) - 1) / 2) * bw
    ax.bar(pos, g["partial_rho"], bw, color=COLORS[name], label=SHOW.get(name, name),
           yerr=[g["partial_rho"] - g["ci_low"], g["ci_high"] - g["partial_rho"]],
           capsize=1.5, error_kw={"lw": 0.7})
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(x)
ax.set_xticklabels(panels)
ax.set_ylabel("partial Spearman rho with tokens/kb")
ax.set_title("Does resolution follow conservation inside each region? (adjusted for GC, repeats, entropy)",
             loc="left", fontweight="bold", fontsize=10)
ax.legend(frameon=False, ncol=5, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
savefig(fig, OUT / "fig_by_region_phylop_only")
log(f"wrote {OUT}")
