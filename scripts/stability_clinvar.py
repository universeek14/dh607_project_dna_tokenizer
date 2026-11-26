"""Analysis 1: how much of the tokenization changes when a ClinVar SNV is applied.

Tokenizers: the 9 from scripts/tokenizer_analysis/common.py. Per variant: bases_retok
(bases covered by alt tokens, as (start, id), absent from ref), n_ref, n_alt.
Input: clinvar_windows.tsv. Output: clinvar_stability_long.tsv. See README.md.
"""
import sys
import time
import pandas as pd

sys.path.insert(0, "/mnt/disk2/vishnu/dna_tokenizers/scripts/tokenizer_analysis")
from common import DNATokenizer

NAMES = ["char", "3mer", "6mer", "NT", "GROVER", "DNABERT-2", "BPE", "EvoLen", "GENA-LM"]

df = pd.read_csv("clinvar_windows.tsv", sep="\t")
df["gc"] = df["ref_window"].apply(lambda s: (s.count("G") + s.count("C")) / len(s))


def pieces(tok, seq):
    toks, starts = tok.encode(seq)
    toks, starts = list(toks), list(starts)
    ends = starts[1:] + [len(seq)]
    return [(s, t, e - s) for s, t, e in zip(starts, toks, ends)]


# quick check against the earlier results (first 300 windows only)
print("CHECK: tokens per kb and mean token length")
for name in NAMES:
    tok = DNATokenizer(name)
    n_tok = total_len = 0
    for seq in df["ref_window"].head(300):
        p = pieces(tok, seq)
        n_tok += len(p)
        total_len += sum(l for s, t, l in p)
    print(f"  {name:10s} tokens/kb {n_tok / 300:7.1f}   mean length {total_len / n_tok:5.2f}")

rows = []
for name in NAMES:
    t0 = time.time()
    tok = DNATokenizer(name)
    for i, r in enumerate(df.itertuples(index=False)):
        a = pieces(tok, r.ref_window)
        b = pieces(tok, r.alt_window)
        a_keys = {(s, t) for s, t, l in a}
        changed = sum(l for s, t, l in b if (s, t) not in a_keys)
        rows.append({"idx": i, "tokenizer": name, "bases_retok": changed,
                     "n_ref": len(a), "n_alt": len(b)})
    print(f"{name} done in {time.time() - t0:.0f} s", flush=True)

out = pd.DataFrame(rows)
meta = df[["group", "label", "gene", "consequence", "gc"]].reset_index().rename(columns={"index": "idx"})
out = out.merge(meta, on="idx")
out.to_csv("clinvar_stability_long.tsv", sep="\t", index=False)

print("\nMean bases re-tokenized per variant, all 12,000:")
print(out.groupby("tokenizer", sort=False)["bases_retok"].mean().round(2))

