# Step 7: Find *PHAS* loci with Phasis

[Back to the skill map](../SKILL.md) · Previous: [06_mirna_validation_mirscore.md](06_mirna_validation_mirscore.md) · Next: [08_structural_rna.md](08_structural_rna.md)

## Goal

21-nt and 24-nt *PHAS* loci, each with the libraries that support it and its register-supported phasiRNAs, plus a separate, weaker *PHAS*-like class.

Phasis alone does the phasing analysis (ShortStack 4 has none). Published *PHAS* catalogs are not used as references; only loci called in the study are annotated.

## The four runs and the comparison

Run **all** Phasis analyses of a species in **one run directory**, **one after another**, each with its own `--outdir`:

```bash
cd runs/phas/phasis        # contains genome/<genome>.fa (+ .fai) as links
phasis -libs tag/*.tag -libformat T -reference genome/<genome>.fa --reference_id_mode preserve \
  -phase 21 -cores 32 --outdir 21_results_individual
phasis ... -phase 24 --outdir 24_results_individual
phasis ... -phase 21 --pool_libraries --outdir 21_results_pooled
phasis ... -phase 24 --pool_libraries --outdir 24_results_pooled
# final set per phase (default --prefer auto keeps the per-library record for shared loci)
phasis-compare 21_results_pooled 21_results_individual --label-a pooled --label-b individual --outdir 21_results_combined
```

- **Why one directory:** Phasis reuses its genome index, processed libraries and `phasis.mem` cache inside a run directory, and the index is expensive to build (25 min for tomato, mostly compression). Two Phasis runs at the same time in one directory are unsafe; a second directory repeats the index cost.
- **Exact matching:** keep `-mismat 0` (the default) and record `-maxhits`.
- `--reference_id_mode preserve` keeps the genome's own sequence names in every output.
- **Individual vs pooled:** individual runs call loci per library, which is replicate-like evidence. Pooled runs merge all libraries, which finds weak loci, but for 24 nt it also picks up dense hc-siRNA clusters with only short phased stretches. Rice pooled 24 nt added 3,756 loci (median 19 kb; the combined set spans about a quarter of the genome).
- **Do not filter pooled 24-nt calls before looking.** First show the strict and exact Howell score distributions for pooled-only loci versus loci also found per library, and let the project lead choose the cut-off. In rice, strict Howell separated the two groups well (AUC 0.90).

## Run time and where it goes

| Run | Tomato (800 Mb, 100 libraries) |
|---|---|
| 21 nt individual, from scratch | 4.3 h at 32 cores |
| 24 nt individual | more than 28 h at 64 cores |
| 24 nt pooled (rice, 374 Mb) | 28.5 h |

Measured slow stages in the 21-nt run: cluster scoring (84 of 257 min), feature assembly (43 min), compressing processed libraries (21 min), plots (21 min), reference compression (20 min), parsing alignments (16 min). The 24-nt run is far slower because about 100× more clusters pass the pre-filter:
- cluster scoring scores every occupied position with three statistical tests, library by library;
- feature assembly then ran for more than 10 hours.

These are tool-level issues: report them to the Phasis developers, do not patch the tool inside an annotation run. Plan the queue so other steps run meanwhile ([04_compute_and_queue.md](04_compute_and_queue.md)).

## Outputs (per `--outdir`)

| File | What it is |
|---|---|
| `<P>_PHAS.gff` | one feature per *PHAS* locus |
| `<P>_calls.tsv` | one row per locus: coordinates, best library (`alib`), Phasis score, Howell scores |
| `<P>_classification_evidence.tsv` | every scored cluster with its final class (PHAS / PHAS-like / non-PHAS) and the reason |
| `<P>_phasiRNAs.tsv` | register-supported phasiRNAs per locus |
| `<P>_all_clusters.tsv` | every candidate cluster with all features |
| `<P>_PHAS_like/` | the same tables for the *PHAS*-like class |
| `<P>_PHAS_locus_plots/`, `*.pdf` | per-locus plots, abundance and Howell score summaries |

How to read them, with tomato numbers: [09_outputs_and_interpretation.md](09_outputs_and_interpretation.md).

## Counting *PHAS* loci

One GFF entry can contain several phased window units (a main unit and secondary units, often neighbouring loci). Report both GFF loci and phased sub-loci, and say which one you count. Rice 21-*PHAS*: 1,011 GFF loci against 1,535-1,667 sub-loci.

## Interpretation limits

- A genomic interval rich in 21- or 24-nt reads is not a *PHAS* locus by that fact alone.
- Reads inside a *PHAS* interval are not automatically phasiRNAs; only register-supported products are.
- Phasing evidence says nothing about mobility, targets or function.
