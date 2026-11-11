# Step 8: Annotate structural RNAs genome-wide

[Back to the skill map](../SKILL.md) · Previous: [07_phas_phasis.md](07_phas_phasis.md) · Next: [09_outputs_and_interpretation.md](09_outputs_and_interpretation.md)

## Goal

Genome-wide tRNA, rRNA, snRNA, snoRNA and organelle-type ncRNA annotations, predicted from the genome first and compared with RNAcentral afterwards.

**Why predict first:** a database-first workflow cannot find what the database lacks or excludes. Maize B73v5 v1 started from RNAcentral and excluded organelle-labelled records, so it contains no nuclear insertions of organelle rRNA. A genome-wide scan of rice found 75 bacterial-type rRNA hits on nuclear chromosomes.

## Run (`scripts/structural/`)

1. Build the Rfam structural subset once per Rfam release: `make_rfam_structural_subset.py` (rRNA, spliceosomal snRNA, C/D and H/ACA snoRNA, scaRNA; 827 families in Rfam 15.1) plus its clan file.
2. `run_structural_scan.sh <genome> <rfam_dir> <outdir> <name> <threads>` runs:
   - **tRNAscan-SE** in eukaryotic mode (`-E`) with EukHighConfidenceFilter. Failing predictions are flagged, not removed.
   - **Infernal `cmscan`** genome-wide: `--cut_ga --rfam --nohmmonly --fmt 2 --clanin`, with Z = 2 × genome size in Mb. Overlapping hits from the same clan compete and only the best is kept.

Tomato took 6.8 h at 32 cores; rice took about 45 min.

## Classification rules

| Case | Class |
|---|---|
| Rfam hit covering ≥ 80% of its model, nuclear-type best model | canonical; RNAcentral support and model coverage as attributes |
| Hit covering < 80% of its model | exploratory (`partial_rfam_hit`); fragments are often repeat-derived |
| Best model a bacterial or archaeal LSU/SSU rRNA | separate class: organelle-type ncRNA (putative NUPT/NUMT if no organelle contigs) |
| RNAcentral placement with no Rfam hit | kept, flagged `rfam_genome_scan=no_hit` for review |

All classes are kept; none is thrown away.

## tRNA specifics

Keep two representations:
- **Genomic:** coordinates, introns included, for the GFF3.
- **Mature:** for read matching. Remove introns from canonically spliced tRNAs, add a terminal CCA only when the genome does not already encode it, and collapse identical mature sequences. Hold low-score intron calls for review.

## Things to review

- Very high copy numbers of one family (rice: 157 snoR71 hits) may be repeat inflation.
- tRNAs with mismatched isotype calls (rice: 108).
- RNAcentral comparison by exact full-length matching is too strict across cultivars; near-identical alignment is a planned improvement.
