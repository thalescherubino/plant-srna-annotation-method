# Step 1: Prepare the reference genome

[Back to the skill map](../SKILL.md) · Next: [02_library_sourcing.md](02_library_sourcing.md)

## Goal

One working copy of the genome FASTA that every tool uses, with its identity and checksums recorded.

## Do

1. **Never edit the original file.** Decompress to a working copy; record SHA-256 checksums of both the original (for example `.fa.gz`) and the copy.
2. **Simplify headers to the first token** (`>1 description...` → `>1`). Tools disagree on how they cut header lines; first-token names keep coordinates identical across all outputs. Keep a table of old → new names if anything changed.
3. **Index it:** `samtools faidx` (needed by miRador, ShortStack, the miRScore input builder) and a Bowtie 1 index (needed by the exact-mapping QC and SRA trimming evaluation). Phasis and ShortStack build their own indexes.
4. **Describe it** in the species record: number of sequences, total length, unplaced contigs, amount of `N`, upper/lower case (soft-masking), and whether chloroplast or mitochondrial contigs are present.
5. **Identify the assembly** from evidence (exact chromosome lengths, naming), not from the file name. Example: a rice genome FASTA supplied without release notes was identified as IRGSP-1.0 in MSU7 packaging from the exact length of Chr1 (43,270,923 bp).

## Why it matters

- Every annotation coordinate refers to this file. A different assembly version, even one with the same chromosome names, silently shifts coordinates.
- **Organelle contigs are usually absent.** Then reads from chloroplast and mitochondrial RNAs (abundant fragments, often > 24 nt) cannot map, and bacterial-like rRNA/tRNA hits on nuclear chromosomes are nuclear insertions (NUPT/NUMT), not organelle genes. Never assume organelles are present.
- **Polyploids** (for example *Coffea arabica*): split the subgenomes and treat each as its own genome. Otherwise reads from homeologous copies multimap and exact-match evidence collapses.

## Check before moving on

- [ ] Checksums of original and working copy recorded
- [ ] Header names are first-token and identical in the FASTA, the `.fai` and the Bowtie index
- [ ] Organelle scope written down (present / absent)
- [ ] Ploidy handled
