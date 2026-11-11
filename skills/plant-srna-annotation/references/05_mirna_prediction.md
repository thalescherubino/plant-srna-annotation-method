# Step 5: Predict miRNA genes (ShortStack 4 and miRador)

[Back to the skill map](../SKILL.md) · Previous: [04_compute_and_queue.md](04_compute_and_queue.md) · Next: [06_mirna_validation_mirscore.md](06_mirna_validation_mirscore.md)

## Goal

Independent lists of candidate *MIR* loci (hairpin, mature, star) from two predictors that use different algorithms, so their agreements and differences can be judged.

## Why two predictors

ShortStack finds sRNA clusters first and then asks whether a cluster looks like a *MIR* locus; miRador starts from inverted repeats in the genome and asks whether reads support a duplex. Each misses some loci the other finds. Keep each predictor's provenance; reconcile afterwards.

## ShortStack 4 (two runs)

```bash
# de novo
ShortStack --genomefile genome.fa --readfile shortstack_condensed/*.fa.gz \
  --dn_mirna --threads 32 --outdir results
# guided by known miRNAs (species miRBase matures, DNA alphabet: U -> T)
ShortStack --genomefile genome.fa --readfile shortstack_condensed/*.fa.gz \
  --dn_mirna --known_miRNAs known_miRNAs.<sp>.dna.fa --threads 32 --outdir results
```

- Use ShortStack **4.1.0 or newer** with the condensed FASTA made from the frozen tag counts.
- `--dn_mirna` turns on de novo *MIR* discovery. **Never add `--nohp`**: it switches miRNA searching off.
- ShortStack 4 has **no phasing analysis**. Its other clusters are reported (for example their overlap with Phasis *PHAS*-like loci) but are not used in later steps.
- Main outputs: `Results.txt` / `Results.gff3` (one row per cluster, `MIRNA` = Y for *MIR* calls), `mir.fasta`, per-library BAM files, `Counts.txt`. See [09_outputs_and_interpretation.md](09_outputs_and_interpretation.md).

Tomato: 107 *MIR* clusters de novo and 120 with known miRNAs, out of about 51,800 clusters.

## miRador

- Configure `miRador.ini`: genome, `libFolder` = the frozen `mirador_libs/` folder, organism code (three letters, for example `sly`), threads, output folder.
- Use a **new** output folder; existing output is overwritten.
- Run unsplit on all libraries: read support is judged across libraries, so splitting the set can change which candidates pass.
- `samtools faidx` is called from `PATH`, so keep the environment active (`mamba run` does this).
- This project uses a Python 3.12 port of miRador, validated to reproduce the original's test outputs byte for byte. Record the code version used.
- Main outputs: `finalAnnotatedCandidates.csv` (candidates with per-library abundances, annotated against miRBase), `preAnnotatedCandidates.csv`, RNAfold images, per-library candidate files.

## Checks before scoring

- [ ] Both ShortStack runs and miRador finished, with their logs and versions saved
- [ ] miRador read the frozen library folder (not an older selection)
- [ ] Candidate counts are plausible (about 100 for ShortStack and a few hundred for miRador in tomato and rice). Zero, or tens of thousands, means a configuration error.
