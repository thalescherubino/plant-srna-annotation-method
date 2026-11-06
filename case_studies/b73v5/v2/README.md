# B73v5 curated sRNA reference: Version 2 plan (structural RNAs)

Decided 2026-11-06. This applies to B73v5 v2 and is the structural-RNA policy for every species in this project. v1 (`B73v5_curated_sRNA_reference_v1_2026-08-27/`) stays frozen; v2 is built in a new versioned directory.

## Why

The v1 structural branch started from RNAcentral records with exact full-length genome mapping and excluded organelle-labelled records. It never searched the genome itself, so it could not find biologically relevant features such as organelle-derived rRNA inserted into nuclear chromosomes (NUPTs/NUMTs). The rice genome-wide scan found 75 bacterial-model rRNA hits (77 including archaeal-model hits) on assembled chromosomes; B73v5 v1 has none, by construction.

## Evidence (same tools and settings for every species)

1. tRNAscan-SE 2.0.13 `-E -H --detail`, then EukHighConfidenceFilter.
2. Genome-wide Infernal 1.1.4 `cmscan --cut_ga --rfam --nohmmonly --fmt 2 --clanin`.
   - Models: the Rfam 15.1 structural subset (827 families: rRNA, snRNA, C/D and H/ACA snoRNA, scaRNA).
   - `-Z` = 2 × genome size in Mb.
   - Non-overlapping hits after clan competition (`olp` `*` or `^`).
3. Comparison with RNAcentral afterwards: exact full-length placements of typed records (as in v1), plus the RNAcentral release's own coordinate file where available.

## Rules

| Case | v2 treatment |
|---|---|
| Rfam hit covering ≥ 80% of its model, nuclear-type best model | canonical; attributes `rnacentral_support=yes/no`, `model_coverage` |
| Rfam hit covering < 80% of its model | exploratory (`partial_rfam_hit`) |
| Rfam hit whose best model is bacterial/archaeal LSU/SSU | separate class `organelle_type_rRNA` (putative NUPT/NUMT); not canonical nuclear rRNA; not in quantification FASTAs |
| v1 RNAcentral placement with no overlapping Rfam hit | kept canonical as in v1, flagged `rfam_genome_scan=no_hit` |
| v1 placement overlapping an Rfam hit | single feature; both evidence sources recorded |
| tRNAs | tRNAscan genome-wide (same as v1) plus high-confidence flags |
| tRNA spliced mature patch (2026-11-05) | included: 37 features spliced (5 mature sequences); 13 low-score features held as in v1 |

## Outputs

- v2 GFF3 (combined and per class).
- Updated structural FASTAs and sequence-equivalence/membership tables.
- A v1→v2 change ledger: added, reclassified and flagged features.
- Checksums, plus validation with the strict checker from 2026-11-05.

miRNA and PHAS components are carried over unchanged from v1.
