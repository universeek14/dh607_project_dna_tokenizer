"""Sample 3,000 variants per (missense/intron x benign/pathogenic) and cut 1-kb windows.

Default: window = 1,000 bp with the variant at index 500 -> clinvar_windows.tsv (stability tests).
--window 1024 --half 512 --out clinvar_windows_1024.tsv gives the 1,024-bp windows (variant at
index 512, the pretraining window length) used for zero-shot scoring. The variant sample is
the same for any window size (sampling happens before cutting, seed 42).
Alt window differs at the variant base only. Windows crossing a chromosome end or containing N
are dropped. Input: clinvar_snv_checked.tsv, hg38.fa.
"""
import argparse

import pandas as pd
from pyfaidx import Fasta

ap = argparse.ArgumentParser()
ap.add_argument("--window", type=int, default=1000, help="total window length")
ap.add_argument("--half", type=int, default=500, help="letters left of the variant (= its index)")
ap.add_argument("--out", default="clinvar_windows.tsv")
args = ap.parse_args()

FASTA = "hg38.fa"
N_PER_GROUP = 3000
HALF = args.half
WINDOW = args.window

df = pd.read_csv("clinvar_snv_checked.tsv", sep="\t")
print("variants in:", len(df))

df = df[df["consequence"].isin(["missense_variant", "intron_variant"])].copy()
df["group"] = df["consequence"].map({"missense_variant": "missense", "intron_variant": "intron"})
print("available per group:")
print(df.groupby(["group", "label"]).size())

parts = []
for (g, lab), grp in df.groupby(["group", "label"]):
    parts.append(grp.sample(n=min(N_PER_GROUP, len(grp)), random_state=42))
sample = pd.concat(parts)
print("sampled per group:")
print(sample.groupby(["group", "label"]).size())

fa = Fasta(FASTA)
rows = []
dropped_edge = dropped_n = 0
for r in sample.itertuples(index=False):
    chrom, pos = r.chrom, int(r.pos)
    start = pos - 1 - HALF
    end = start + WINDOW
    if start < 0 or end > len(fa[chrom]):
        dropped_edge += 1
        continue
    ref_w = fa[chrom][start:end].seq.upper()
    if "N" in ref_w:
        dropped_n += 1
        continue
    assert ref_w[HALF] == r.ref
    alt_w = ref_w[:HALF] + r.alt + ref_w[HALF + 1:]
    rows.append({"chrom": chrom, "pos": pos, "ref": r.ref, "alt": r.alt,
                 "label": r.label, "group": r.group, "consequence": r.consequence,
                 "gene": r.gene, "ref_window": ref_w, "alt_window": alt_w})

out = pd.DataFrame(rows)
out.to_csv(args.out, sep="\t", index=False)
print("dropped near chromosome ends:", dropped_edge)
print("dropped for unknown letters (N):", dropped_n)
print("final windows:", len(out))
print(out.groupby(["group", "label"]).size())
