# Step 9: Read and interpret the outputs

[Back to the skill map](../SKILL.md) · Previous: [08_structural_rna.md](08_structural_rna.md) · Next: [10_release.md](10_release.md)

For each output: what one row is, the columns that matter, typical values, and what the result does **not** show. Numbers are preliminary rice and tomato values for orientation, not targets.

## Library QC tables (step 3)

- **One row = one library.** Key columns: reads kept, length percentages (`pct21`, `pct24`), `5U_21nt`, `5A_24nt`, exact-mapping rate, top-sequence share, flags.
- **A good plant library:** peaks at 21 and/or 24 nt, 5′U ≥ 0.4 at 21 nt, 5′A ≥ 0.4 at 24 nt, most 18-30-nt reads mapping exactly.
- **Doesn't show:** tissue identity. Metadata can be wrong, and an odd profile can be real biology (pollen and seeds differ strongly from leaves).

## ShortStack `Results.txt`

- **One row = one sRNA cluster.** `Locus`, `Reads`, `DistinctSequences`, `FracTop` (fraction of reads on the top strand), `Strand`, `MajorRNA` (most abundant sequence), counts per length (`21`…`24`), `DicerCall` (dominant size, or `N` when there is none), `MIRNA` (Y = called a *MIR* locus). The known-miRNA run adds `known_miRNAs`.
- **Typical:** tomato about 51,800 clusters; 107 `MIRNA = Y` de novo, 120 with known miRNAs. Most clusters are 24-nt-dominant.
- **Doesn't show:** a `DicerCall` of 24 is not an hc-siRNA label, and `N` clusters may be degradation fragments. Both are reported, not used further, until reviewed.

## miRador `finalAnnotatedCandidates.csv`

- **One row = one candidate mature miRNA** (both arms of a duplex appear as separate rows). Columns: chromosome, strand, positions and sequences of miR and star, number of genome hits, then abundance per library, abundance of 1-nt variants and the precursor's total.
- Candidates are named `candidate-<n>_<k>-5p/3p`, or by their miRBase family when they match one.

## miRScore `<set>_miRScore_results.csv`

- **One row = one hairpin.** Columns: `mSeq`/`msSeq` (mature/star), positions and lengths, `result` (Pass/Fail), `flags`, `mismatches` and `bulges` in the duplex, `mReads`/`msReads`/`totReads`, `numLibraries`, `precision` (%).
- **Flags seen in tomato miRador candidates (481 hairpins):**

| Flag | Meaning | Pass / Fail |
|---|---|---|
| `NA` | no problem found | 273 / 0 |
| `23/24 nt miRNA` (+ `star`) | mature or star is 23-24 nt: unusual for miRNAs; can be a real 24-nt miRNA or a hairpin-derived siRNA | 69 / 117 |
| `Precision less than 75%` | reads are scattered over the hairpin, not concentrated on the duplex | 0 / 121 |
| `No 2nt 3' overhang` | duplex lacks the Dicer signature | 0 / 13 |
| `More than 5 mismatches in duplex` | pairing too poor | 0 / 4 |
| `No mature or star reads detected` | no read support for one strand | 0 / 4 |
| `Less than 10 reads in a single library` | support never reaches the minimum in any library | 0 / 2 |
| `Hairpin structure invalid` | the precursor does not fold into a usable stem | 0 / 1 |

- **Passing tomato miRador hairpins:** median precision 94%.
- `alt_mirna_results.csv` lists rescued alternative duplexes (45 for tomato miRador, all Pass). A rescued Pass means the locus is supported but the annotated duplex is not the dominant one; review before replacing the original annotation.
- **Doesn't show:** experimental validation, targets or function. "Pass" means "meets plant miRNA annotation criteria on these libraries". A Fail can also mean the tissue that expresses that miRNA is missing from the library set.

## Phasis tables (per phase, per run)

- **`<P>_calls.tsv`, one row = one *PHAS* locus in one supporting library:** `identifier` (`chr:start..end`, the locus), `alib` (the library), `phasis_score` (capped at 300), `Peak_Howell_score`, `Peak_Howell_score_strict`, `Howell_exact_support_score`, `Howell_origin_class`. **Count distinct `identifier` values for loci**; the number of rows per locus is how many libraries support it.
  - Tomato 21-nt individual: 1,006 rows = 87 loci (the GFF has 87 lines), supported by 96 libraries in total. 24-nt individual: 898 rows = 366 loci in 25 libraries. 21-nt pooled: 17 loci, broader (median 24 kb against 4 kb per library), because pooled reads chain neighbouring clusters together.
  - The `Howell_origin_class` counts below are per row, not per locus.
  - `Howell_origin_class`: `coherent_extension` 769, `mixed_extension_and_ambiguity` 178, `unique_origin` 53, `ambiguous_origin` 6. This describes whether the phased signal grows coherently from one origin or competes with alternative registers.
- **`<P>_classification_evidence.tsv`, one row = one scored cluster per library**, with `final_class` and `evidence_reason`.
  - Tomato 21 nt (rows, i.e. cluster × library): 35,691 non-*PHAS*, 1,006 *PHAS*, 622 *PHAS*-like.
  - *PHAS*-like reasons: `insufficient_exact_support`, `weak_exact_support`, `low_score_crowded_window_context`, `weak_scaffold_alternative_context`.
- **`<P>_phasiRNAs.tsv`, one row = one phasiRNA sequence at a locus:** `strand` (w = top/Watson, c = bottom/Crick), `observed_pos` vs `expected_register_pos`, `register_class` (`core_exact`, `core_offset`, `extended_exact`), `window_unit_role` (`main_hpsp`, `main_partner`, `overlapping_alternative`, `other_local_peak`), `abun`, `tag_seq`, `hits`.
  - Tomato 21 nt: 28,371 phasiRNA rows.
  - Sequences are stored reference-forward, as the read appears on the genome's top strand, also for minus-strand loci. Do not reverse-complement them again.
- **`<P>_PHAS.gff`:** one line per distinct locus (merged across libraries) with summary attributes (`complexity`, `strand_bias`, `ratio_abund_len_phase`).
- **Howell score distributions:** compare pooled-only loci with loci also found per library before choosing any cut-off ([07_phas_phasis.md](07_phas_phasis.md)).
- **Doesn't show:** a *PHAS* locus is not necessarily a protein-coding gene or a known *TAS* locus, and phasiRNAs are not shown to move or act.

## Structural scan

- **tRNAscan-SE:** one row per tRNA, with isotype, anticodon, intron, score, pseudogene call and high-confidence flag.
- **cmscan table:** one row per hit, with family, clan, score, E-value and model coordinates (coverage is computed from these).
- **Rice:** 733 tRNAs (538 high confidence); 1,046 Rfam hits (562 C/D snoRNA, 336 rRNA, 83 snRNA, 65 H/ACA snoRNA); 75 organelle-type rRNA hits.
- **Doesn't show:** expression. A hit is a gene model, not evidence that the gene is transcribed.

## Cross-class review

The same sequence can come from several loci or classes (for example a tRNA fragment that also lies in a *PHAS* interval). Keep every membership and do not count abundance twice. Decide precedence rules explicitly and write them down.
