"""Heatmaps of the by-region partial correlations (all features; conservation only).

Input/output folder: B_by_region/ (argv[1]).
"""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = sys.argv[1] if len(sys.argv) > 1 else "B_by_region"
res = pd.read_csv(f"{folder}/partial_correlations_by_region.csv")
res["sig"] = (res.ci_low > 0) | (res.ci_high < 0)

TOK = ["GROVER", "DNABERT-2", "BPE", "EvoLen", "GENA-LM"]
REG = ["All regions", "promoter", "exon", "CTCF", "enhancer", "intron", "random", "repeat"]
FEAT = [("phylop_mean", "Conservation (mean phyloP)"), ("gc", "GC fraction"),
        ("rep_frac", "Repeat fraction"), ("kmer3_entropy", "3-mer entropy")]
SHOW = {"BPE": "Standard BPE"}


def grid(feature):
    s = res[res.feature == feature]
    val = s.pivot(index="region", columns="tokenizer", values="partial_rho").loc[REG, TOK]
    sig = s.pivot(index="region", columns="tokenizer", values="sig").loc[REG, TOK]
    return val, sig


def heat(ax, feature, title, lim):
    val, sig = grid(feature)
    im = ax.imshow(val.values, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(TOK)))
    ax.set_xticklabels([SHOW.get(t, t) for t in TOK], rotation=30, ha="right", fontsize=9)
    ax.set_yticks(range(len(REG)))
    ax.set_yticklabels(REG, fontsize=9)
    ax.axhline(0.5, color="black", lw=1.2)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    for i in range(len(REG)):
        for j in range(len(TOK)):
            v = val.values[i, j]
            txt = f"{v:.2f}" + ("*" if sig.values[i, j] else "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 0.6 * lim else "black")
    return im


NOTE = ("* = 95% range (resampled over 10-Mb genome blocks within each region) excludes zero. "
        "Top row is all regions pooled.")

# Figure 1: all four features on one shared color scale
fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
for ax, (f, t) in zip(axes.flat, FEAT):
    im = heat(ax, f, t, 0.9)
fig.colorbar(im, ax=axes, shrink=0.6, label="partial Spearman rho with tokens/kb")
fig.suptitle("What sets tokenization resolution? Each feature adjusted for the other three, by region",
             fontweight="bold")
fig.text(0.5, -0.01, NOTE + " Same color scale in all four panels.", ha="center", fontsize=8, color="gray")
fig.savefig(f"{folder}/heatmap_all_features.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{folder}/heatmap_all_features.pdf", bbox_inches="tight")

# Figure 2: conservation only, on its own, tighter color scale
fig, ax = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
im = heat(ax, "phylop_mean", "Does resolution follow conservation? (adjusted for GC, repeats, entropy)", 0.2)
fig.colorbar(im, ax=ax, label="partial Spearman rho with tokens/kb")
fig.text(0.5, -0.02, NOTE + " Color scale is narrower than in the four-panel figure.",
         ha="center", fontsize=8, color="gray")
fig.savefig(f"{folder}/heatmap_conservation.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{folder}/heatmap_conservation.pdf", bbox_inches="tight")
print("saved heatmaps in", folder)
