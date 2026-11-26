# dh607_project_dna_tokenizer

Computational Multiomics (DH607) course project: testing EvoLen, a DNA tokenizer we expected to work well, against eight other tokenizers. This repo holds our tokenizer analysis scripts and results. Matched small-BERT pretraining was done separately.

## Motivation

An adaptive tokenizer could use finer resolution in informative genomic regions and coarser resolution elsewhere. Before building one, we asked how existing tokenizers behave on the genome, and whether EvoLen, which uses evolutionary information, behaves differently. This repo does not build an adaptive tokenizer.

## What we did

1. **Tokenizer characterization.** Nine tokenizers (character, 3-mer, 6-mer, NT-v2, GROVER, DNABERT-2, BPE, EvoLen, GENA-LM) on 28,000 one-kb windows of human hg38 from seven region classes. We measured token density, token lengths, links to GC, repeats and conservation, alignment with exon and splice boundaries, and response to single-base changes.
2. **ClinVar stability test.** For 12,000 real ClinVar variants (3,000 each of missense or intronic, pathogenic or benign), we cut a 1,000-letter hg38 window around the variant, made a copy with the variant letter changed, and counted how many letters each tokenizer re-cuts. Pathogenic and benign variants were compared with GC adjustment and gene-resampled confidence ranges.
3. **Matched pretraining (separate work).** One small BERT per tokenizer (about 19M parameters, 40,000 steps, same data and settings) for nine tokenizers, including a two-by-two grid that separates EvoLen's two ingredients.
4. **Exploratory link.** Tokenizer instability compared with the pretrained models' scores.

## Main observations

- Learned tokenizers follow GC and repeat content, not necessarily conservation, and rarely place cuts at exon or splice boundaries.
- In case of adaptive tokenizers, EvoLen re-cuts the most letters after a one-letter change (about 12 against about 9 for BPE).
- Pathogenic variants are re-cut slightly less than benign ones in every learned tokenizer, and this does not depend on conservation.
- More unstable tokenizers tend to score lower in the pretrained models. This is exploratory.

## Run information

Scripts contain absolute paths from our server and import helper code from the shared project, so edit the paths before running. Large files (genome, raw ClinVar, generated tables) are not included.

## Data sources

ClinVar (NCBI, GRCh38 VCF), UCSC hg38 with phyloP 100-way, ENCODE SCREEN, JASPAR 2024, GENCODE v44. Tokenizers from Hugging Face (GROVER, DNABERT-2, GENA-LM, NT-v2) and the EvoLen repository (arXiv 2604.08698).

## Limitations

- Most analyses measure the tokenizer input, not model predictions.
- ClinVar covers clinically studied variants, mostly inside genes.
- We do not know why pathogenic variants are re-cut less.
- Models are small, with one training run per tokenizer (two for four of them).

## Credits

Analysis scripts: Sanskriti Kadam. Pretraining, downstream evaluation and shared tokenizer code: Vishnu A M (add link). EvoLen belongs to its authors.
