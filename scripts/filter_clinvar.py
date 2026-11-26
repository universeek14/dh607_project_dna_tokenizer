"""Filter the ClinVar GRCh38 VCF to labelled, review-filtered SNVs.

Input: raw/clinvar.vcf.gz. Keeps single-base A/C/G/T SNVs on chr1-22/X with a clear
(likely) benign = 0 or (likely) pathogenic = 1 significance and review status
"criteria provided" or better. Consequence = first MC entry; gene = first GENEINFO symbol.
Output: clinvar_snv_labeled.tsv. Run: python filter_clinvar.py
"""
import gzip

PATH = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
GOOD_REVIEW = ("criteria_provided", "reviewed_by_expert_panel", "practice_guideline")
CHROMS = {str(i) for i in range(1, 23)} | {"X"}

total = snv = labeled = 0
kept = {0: 0, 1: 0}

out = open("clinvar_snv_labeled.tsv", "w")
out.write("chrom\tpos\tref\talt\tlabel\tclnsig\treview\tconsequence\tgene\n")

with gzip.open("raw/clinvar.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        total += 1
        p = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, info = p[0], p[1], p[3], p[4], p[7]
        if chrom not in CHROMS:
            continue
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
            continue
        d = {}
        for item in info.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                d[k] = v
        if d.get("CLNVC") != "single_nucleotide_variant":
            continue
        snv += 1
        sig = d.get("CLNSIG", "")
        if sig in PATH:
            label = 1
        elif sig in BENIGN:
            label = 0
        else:
            continue
        labeled += 1
        rev = d.get("CLNREVSTAT", "")
        if not rev.startswith(GOOD_REVIEW):
            continue
        mc = d.get("MC", "NA").split(",")[0].split("|")[-1]
        gene = d.get("GENEINFO", "NA").split(":")[0]
        out.write("\t".join(["chr" + chrom, pos, ref, alt, str(label), sig, rev, mc, gene]) + "\n")
        kept[label] += 1

out.close()
print("all records:", total)
print("single-letter variants:", snv)
print("clearly labeled:", labeled)
print("kept after review filter - benign:", kept[0], " pathogenic:", kept[1])
