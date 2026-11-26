"""Analysis 3: the Analysis 1 stability test for the 9 pretraining tokenizers.

Tokenizers from pretrain_small_bert/tokenizers/tok_registry.py.
Input: clinvar_windows.tsv. Output: clinvar_stability_long_ablation.tsv.
"""
import sys
import time
import pandas as pd

sys.path.insert(0, "/mnt/disk2/vishnu/dna_tokenizers/pretrain_small_bert/tokenizers")
from tok_registry import get_tokenizer

NAMES = ["char", "kmer6", "bpe", "len2_only", "strat_only", "evolen",
         "evolen_rebuilt", "random_strat", "func_strat"]

df = pd.read_csv("clinvar_windows.tsv", sep="\t")
df["gc"] = df["ref_window"].apply(lambda s: (s.count("G") + s.count("C")) / len(s))


def pieces(tok, seq):
    ids, starts = tok.encode(seq)
    ids, starts = list(ids), list(starts)
    ends = starts[1:] + [len(seq)]
    return [(s, t, e - s) for s, t, e in zip(starts, ids, ends)]


rows = []
for name in NAMES:
    t0 = time.time()
    tok = get_tokenizer(name)
    n_tok = 0
    for i, r in enumerate(df.itertuples(index=False)):
        a = pieces(tok, r.ref_window)
        b = pieces(tok, r.alt_window)
        a_keys = {(s, t) for s, t, l in a}
        changed = sum(l for s, t, l in b if (s, t) not in a_keys)
        n_tok += len(a)
        rows.append({"idx": i, "tokenizer": name, "bases_retok": changed,
                     "n_ref": len(a), "n_alt": len(b)})
    print(f"{name} done in {time.time() - t0:.0f} s, mean tokens/kb {n_tok / len(df):.1f}", flush=True)

out = pd.DataFrame(rows)
meta = df[["group", "label", "gene", "consequence", "gc"]].reset_index().rename(columns={"index": "idx"})
out = out.merge(meta, on="idx")
out.to_csv("clinvar_stability_long_ablation.tsv", sep="\t", index=False)

print("\nMean bases re-tokenized per variant:")
print(out.groupby("tokenizer", sort=False)["bases_retok"].mean().round(2))
print("\nBy label (0 = benign, 1 = pathogenic):")
print(out.groupby(["tokenizer", "label"], sort=False)["bases_retok"].mean().round(2).unstack())
