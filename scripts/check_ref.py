"""Check every ClinVar REF base against the hg38 FASTA; keep matching rows.

Input: clinvar_snv_labeled.tsv and the FASTA path (argv[1]).
Output: clinvar_snv_checked.tsv. Run: python check_ref.py hg38.fa
"""
import sys
from pyfaidx import Fasta

fa = Fasta(sys.argv[1])
print("first sequence names:", list(fa.keys())[:5])

n = bad = 0
with open("clinvar_snv_labeled.tsv") as f, open("clinvar_snv_checked.tsv", "w") as out:
    out.write(f.readline())
    for line in f:
        p = line.rstrip("\n").split("\t")
        chrom, pos, ref = p[0], int(p[1]), p[2]
        base = fa[chrom][pos - 1:pos].seq.upper()
        n += 1
        if base == ref:
            out.write(line)
        else:
            bad += 1

print("checked:", n, " mismatches:", bad)

