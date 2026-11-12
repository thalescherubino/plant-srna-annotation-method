# Step 9: Read and interpret the outputs

[Back to the skill map](../SKILL.md) · Previous: [08_structural_rna.md](08_structural_rna.md) · Next: [10_release.md](10_release.md)

For each output: what one row is, the columns that matter, typical values, and what the result does **not** show. Numbers are preliminary rice and tomato values for orientation, not targets.

## Library QC tables (step 3)

- **One row = one library.** Key columns: reads kept, length percentages (`pct21`, `pct24`), `5U_21nt`, `5A_24nt`, exact-mapping rate, top-sequence share, flags.
- **A good plant library:** peaks at 21 and/or 24 nt, 5′U ≥ 0.4 at 21 nt, 5′A ≥ 0.4 at 24 nt, most 18-30-nt reads mapping exactly.
- **Doesn't show:** tissue identity. Metadata can be wrong, and an odd profile can be real biology (pollen and seeds differ strongly from leaves).

## ShortStack `Results.txt`

- **One row = one sRNA cluster.** `Locus`, `Reads`, `DistinctSequences`, `FracTop` (fraction of reads on the top strand), `Strand`, `MajorRNA` (most abundant sequence), counts per length (`21`…`24`), `DicerCall` (dominant size, or `N` when there is none), `MIRNA` (Y = called a *MIR* locus). The known-miRNA run adds `known_miRNAs`.
- **Typical:** tomato about 51,800 clusters; 107 `MIRNA = Y` in the discovery run, 120 with known miRNAs. Most clusters are 24-nt-dominant.
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

Phasis explains these files in its README ("Understand the Outputs") and in `docs/scoring.md` and `docs/phasirna_output.md`. Column meanings below follow those documents; where they are silent, the Phasis code was read and the source is named.

- **`<P>_calls.tsv`, one row = one *PHAS* locus in one supporting library:** `identifier` (`chr:start..end`, the locus), `achr`/`start`/`end`, `alib` (the library; `ALL_LIBS` in a pooled run), `phasis_score`, `Peak_Howell_score`, `Peak_Howell_score_strict`, `Howell_exact_support_score`, `Howell_origin_class`. **Count distinct `identifier` values for loci**; the number of rows per locus is how many libraries support it.
  - Tomato 21-nt individual: 1,006 rows = 87 loci (the GFF has 87 lines), supported by 96 libraries in total. 24-nt individual: 898 rows = 366 loci in 25 libraries. 21-nt pooled: 17 loci, broader (median 24 kb against 4 kb per library), because pooled reads chain neighbouring clusters together.
- **The scores** (Phasis `docs/scoring.md`). The Phasis score and the Howell scores are different measures, not two names for one score.
  - `phasis_score`: statistical evidence for phasing, -log10 of the cluster's window P values combined by Fisher's method, capped at 300. Values at 300 hit the cap; they are not identical evidence.
  - `Peak_Howell_score` (the main Howell score): support for regularly spaced phase-length reads, from their abundance and the number of occupied phase cycles in a window of ten cycles (210 nt for 21 nt, 240 nt for 24 nt). It tolerates reads 1 nt off the expected position. Phasis's `-min_Howell_score` filter (default 12.5) applies to this score.
  - `Peak_Howell_score_strict`: a separate scan that counts only reads at exact register positions. Its peak may come from a different window or register than the main score.
  - `Howell_exact_support_score`: exact-register support inside the main score's peak window. Not the same as the strict score.
- **`Howell_origin_class`** (not described in the README; from the code, `_summarize_peak_howell_ambiguity` in `phasis/stages/feature_assembly.py`): looks at other windows on the same strand that overlap the peak window and reach at least 90% of its exact support.
  - `unique_origin`: no such window. `coherent_extension`: all such windows share the peak's register (the phased signal extends). `ambiguous_origin`: all are in another register (competing registers). `mixed_extension_and_ambiguity`: both kinds.
  - Tomato 21-nt individual, counted per row, not per locus: `coherent_extension` 769, `mixed_extension_and_ambiguity` 178, `unique_origin` 53, `ambiguous_origin` 6.
- **`<P>_classification_evidence.tsv`, one row = one classified cluster in one library** (`identifier`, `alib`, `cID`), with `final_class` (`PHAS`, `PHAS-like` or `non-PHAS`) and `evidence_reason`. The README calls this the register-resolved evidence; the reason values are not listed there and come from the code (`phasis/stages/classify.py`, where the thresholds are):

| `final_class` | `evidence_reason` | Plain meaning | Tomato 21 nt (rows) |
|---|---|---|---|
| *PHAS* | `pass` | passed the classifier and every exact-support check | 1,006 |
| *PHAS*-like | `weak_exact_support` | some exact-register support, but little | 368 |
| *PHAS*-like | `low_score_crowded_window_context` | modest main Howell score among many crowded competing windows | 244 |
| *PHAS*-like | `weak_scaffold_alternative_context` | medium main Howell score, weak exact or strict support, and competing or ambiguous context | 10 |
| non-PHAS | `classifier_non_phas` | rejected by Phasis's classifier | 35,229 |
| non-PHAS | `insufficient_exact_support` | no exact-register support at all | 462 |

- **The class is decided per library, not per locus.** The same locus can be *PHAS* in one library and *PHAS*-like in another: in tomato 21 nt, 62 of the 99 *PHAS*-like loci are also *PHAS* in some other library. The *PHAS*-like GFF therefore overlaps the *PHAS* GFF; remove those identifiers before reporting "*PHAS*-like only" loci.
- **`<P>_phasiRNAs.tsv`, one row = one read assignment** (a sequence in one phased unit of one locus in one library), not one row per distinct sequence (Phasis `docs/phasirna_output.md`). Only reads whose length equals the phase are included. Columns: `identifier`, `cID`, `alib`, `phase`, `window_unit_id` (`unit_main` or `unit_secondary_<n>`, local to the locus and library), `window_unit_role`, `strand` (`w` forward, `c` reverse), `observed_pos` vs `expected_register_pos`, `register_class`, `abun`, `tag_seq`, `hits`.
  - `register_class`: `core_exact` (read at an expected position of the core register), `core_offset` (no read at the exact position, so a read 1 nt away; the more abundant neighbour if both exist), `extended_exact` (exact position in a supported extension beyond the core window). For an exact-only set, drop `core_offset`.
  - `window_unit_role`: the main unit has `main_hpsp` (the anchor window), `main_partner` (compatible opposite-strand support) and `main_extension`; secondary units are `overlapping_alternative` or `other_local_peak`. `other_local_peak` does not by itself mean a physical gap; check positions and the plot.
  - Tomato 21 nt: 28,371 rows (`extended_exact` 14,242, `core_exact` 9,449, `core_offset` 4,680; no `main_extension` rows).
  - **Do not sum `abun` over rows:** one read can be assigned to several units, positions or libraries.
  - **Orientation:** `tag_seq` is written as aligned, i.e. reference-forward (as on the genome's forward strand). `w` rows are the sRNA as written. `c` rows are the reverse complement of the sRNA; reverse-complement them once to get the sRNA 5′→3′, for example for target prediction (Phasis `docs/phasirna_output.md`, "Preparing sequences for target prediction"). Never flip `w` rows, and never flip a sequence twice. State the orientation of any released phasiRNA FASTA (the B73 v1 release kept the reference-forward form).
- **`<P>_PHAS.gff`:** one line per distinct `identifier`, on `achr` from `start` to `end`. The score column is `phasis_score` and the strand is `.`. The attributes (`id`, `complexity`, `strand_bias`, `log_clust_len_norm_counts`, `ratio_abund_len_phase`) come from the first library row of that locus, not from a merge across libraries (code: `write_gff` in `phasis/stages/output.py`).
- **Howell score distributions:** compare pooled-only loci with loci also found per library before choosing any cut-off ([07_phas_phasis.md](07_phas_phasis.md)).
- **Doesn't show:** a *PHAS* locus is not necessarily a protein-coding gene or a known *TAS* locus, and phasiRNAs are not shown to move or act.

## Structural scan

- **tRNAscan-SE:** one row per tRNA, with isotype, anticodon, intron, score, pseudogene call and high-confidence flag.
- **cmscan table:** one row per hit, with family, clan, score, E-value and model coordinates (coverage is computed from these).
- **Rice:** 733 tRNAs (538 high confidence); 1,046 Rfam hits (562 C/D snoRNA, 336 rRNA, 83 snRNA, 65 H/ACA snoRNA); 75 organelle-type rRNA hits.
- **Doesn't show:** expression. A hit is a gene model, not evidence that the gene is transcribed.

## Cross-class review

The same sequence can come from several loci or classes (for example a tRNA fragment that also lies in a *PHAS* interval). Keep every membership and do not count abundance twice. Decide precedence rules explicitly and write them down.
