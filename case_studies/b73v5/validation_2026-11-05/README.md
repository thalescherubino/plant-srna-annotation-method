# B73v5 curated sRNA reference v1 — GFF3/FASTA validation against the local genome

> Copied into this repository from the 2026-11-05 run. The scripts in this folder are the ones that produced it; `run_validation.sh` records the server paths used.

Run: 2026-11-05 · Directory: `/workdir/thc47/gff_validation_2026-11-05/` · Reproduce: `scripts/run_validation.sh`

Validation code here is **new** (`scripts/`). It is not the original annotation-building code, which is not in the release
(`provenance/annotation_build_manifest.tsv` names `scripts/build_b73v5_unified_srna_gff3.py`, located on the original HPC).
`ms28_sharing_package_final_2026-09-14/code/` holds downstream DA/figure code (`phasis_srna_da`), not reference construction code.
All inputs were read-only; the genome `.fai` was written beside a symlink in `genome/`.

## Inputs

| Input | SHA-256 | Note |
|---|---|---|
| `annotation/B73v5_sRNA_reference_v1.combined.gff3` | `f89c85ab…6766bf` | matches release `SHA256SUMS` and build manifest |
| 12 release FASTAs, 6 metadata tables | see `results/strict/input_checksums.tsv` | all 18 match release `SHA256SUMS` |
| `compare_maize/Zm-B73-REFERENCE-NAM-5.0.fa` | `52f06632…de2c11` | **differs** from the report's `a40ee358…ded2c2` (see caveats) |
| Bowtie1 SAMs (21/24 phasiRNAs vs. this genome, 2026-11-04) | in `input_checksums.tsv` | used read-only for placement counts |

Software: Python 3.9.25, samtools 1.20, bedtools 2.29.2; Phasis source inspected at `v2.8.4-2-g90a3a0c`.

## What was tested

For each of the 87,606 combined-GFF3 features (each gets exactly one row in `results/strict/features.all_classes.tsv`):

1. **Coordinate compatibility**: interval lies on a contig of this genome, within bounds.
2. **Feature-specific sequence identity**: the feature's own sequence is resolved through release identifiers, then compared with the genome at that interval:
   - phasiRNA `exon`: GFF `Name` (identical to the sequence embedded in `ID`).
   - `miRNA`: `reference_product_id` → `miRNA_bona_fide.mature_and_star.fa`.
   - MIR `gene`: same ID in `precursors.fa`, compared at the representative-model coordinates.
   - rRNA/snRNA/snoRNA/tRNA: `source_feature_id` → `structural_ncRNA.sequence_membership.tsv` → FASTA `reference_id`.
3. **Links**: FASTA record, `sequence_equivalence.tsv` membership, and parent containment. Also:
   - phasiRNAs: the Phasis observation rows (locus identifier, w/c strand, position, tag; row/library/abundance totals must equal the GFF attributes).
   - miRNA and structural RNAs: the coordinate must be among the release's own listed exact placements.
   - phasiRNAs: the coordinate must be among the Bowtie placements of the feature's FASTA record.

## Results (canonical tier)

| Class | Features | Coord. valid | Exact for own ID | RC only | Other ID | No match | Not evaluable | Multi-placement |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 21-PHAS phasiRNA (`exon`) | 30,792 | 30,792 | **30,792** | 0 | 0 | 0 | 0 | 5,028 |
| 24-PHAS phasiRNA (`exon`) | 27,331 | 27,331 | **27,331** | 0 | 0 | 0 | 0 | 9,242 |
| MIR locus (`gene`) | 220 | 220 | **220** | 0 | 0 | 0 | 0 | – |
| MIR precursor model | 465 | 465 | 220 (representative) | 0 | 0 | 0 | 245¹ | – |
| mature/star (`miRNA`) | 930 | 930 | **930** | 0 | 0 | 0 | 0 | 4 |
| tRNA | 2,088 | 2,088 | **2,088** | 0 | 0 | 0 | 0 | – |
| rRNA placements | 16,875 | 16,875 | **16,875** | 0 | 0 | 0 | 0 | 15,869 |
| snRNA | 499 | 499 | **499** | 0 | 0 | 0 | 0 | 33 |
| snoRNA | 780 | 780 | **780** | 0 | 0 | 0 | 0 | 148 |
| PHAS loci (Phasis) | 426 + 437 | all | – | – | – | – | locus | – |
| PHAS loci, locus-only | 13 + 15 | all | – | – | – | – | locus | – |
| phasing windows; tRNA exons | 3,925; 2,138 | all | – | – | – | – | child | – |

¹ No per-model precursor sequence exists in the release for non-representative models.

All link checks passed for every evaluated feature (0 failures in `summary.per_class.tsv`):

- All 58,123 phasiRNA children trace to Phasis observation rows, and their `source_product_rows`, `source_library_count` and `source_abundance_sum` agree exactly.
- 28 locus-only PHAS additions (13 × 21-nt, 15 × 24-nt) have valid coordinates and zero children. No children were inferred.
- 20 PHAS loci have GFF spans wider than their curated spans; all 20 are declared via `gff_span_expanded_for_phasiRNA_children=true`.
- Every record of every canonical FASTA is referenced by ≥1 GFF feature (`fasta_coverage.tsv`).

**Negative control.** Shifting every phasiRNA and mature/star interval by ±1 nt drops own-ID exact matches from 58,123 → 0 and 930 → 0, so the check is not vacuous.

## phasiRNA orientation convention: resolved

- **Data.** All 58,123 phasiRNA exons match their own `Name` in **reference-forward** orientation, including all 28,294 minus-strand exons (`phasiRNA_orientation_crosstab.tsv`).
- **Code.** In Phasis at the release commit `00fb8c0` (present in the local git history), `phasis/stages/sam_parsing.py` sets `tag_seq` from the SAM `SEQ` field verbatim and sets the strand from flag 16. No reverse-complement step exists anywhere in `phasis/` at that commit, and `sam_parsing.py` is unchanged through v2.8.4. SAM stores reverse-strand reads reference-forward, so Phasis tags for `c`-strand phasiRNAs are reference-forward.

So the `21_PHAS`/`24_PHAS` FASTA and GFF `Name` hold the **genomic plus-strand sequence**. For `-` strand features, the biological (read) sequence is the reverse complement. This explains the earlier stranded bedtools result (15,381 same-strand, only 101 of them minus-strand). Those 101 minus-strand matches are features whose reverse complement happens to be a separate phasiRNA record (`BOTH_ORIENTATIONS` in the class-level helper). Do not edit strands. Downstream read counting must account for this, e.g. by mapping reads against the reverse complement or matching either orientation.

## Findings that need a decision (not errors in coordinates)

1. **Intron-containing tRNAs (50 features, 39 reference sequences).** `tRNA.genomic.fa` holds the **full genomic span including the intron**, and v1 `tRNA.CCA_aware_quantification.fa` inherited it. That is fine for coordinates but wrong for mature-tRNA quantification.
   - **Decision (2026-11-05):** keep the intron in the GFF3 and the genomic FASTA, and remove it from the mature CCA-aware sequence.
   - **Implemented** as a non-destructive patch: `/workdir/thc47/B73v5_curated_sRNA_reference_v1_patch_tRNA_spliced_2026-11-05/`. See its README.
   - All 50 introns are at the canonical position (anticodon ends 1 nt before the intron).
   - 35 references (37 features) were spliced and collapse to 5 mature sequences (881 → 851 quantification records). An explicit-intron sidecar GFF3 is provided.
   - 4 references (13 features) are held at v1 pending review: low tRNAscan scores (30–37), including 11 scaffold Leu-TAA copies whose "spliced" 70-nt form lacks the Leu variable arm.
2. **Exploratory unsupported-MIR products (250).** Results by product:
   - 179 lie exactly within their own coordinate model.
   - 33 have no coordinate model.
   - **38 lie within their model only with 1 mismatch (33) or 2 mismatches (5).** All 38 come from non-exact projections (minimap2, bowtie_v3 remediation, or local alignment). They keep the source miRBase/S4 sequence rather than the B73 allele, unlike canonical models, which were replaced with the oriented genomic allele.

   PHAS-like candidate products: all 177 (21-nt) and 2,315 (24-nt) lie inside their own loci, reference-forward.

   **Open research question (deferred, 2026-11-05).** Mature-miRNA mismatch tolerance should probably be allowed only *across* species, not within one. This needs evidence, measured in three separate ways:
   - **(a) Source-vs-genome allele at the same locus.** This is the case for these 38: a projection or accession artefact, not biology.
   - **(b) Within-species paralog variation among bona fide MIR family members.** Measurable from the 220 canonical loci in `miRNA_bona_fide.products.tsv`.
   - **(c) Cross-species divergence.**

   The exploratory candidates may be non-miRNAs or young, evolving MIR genes. Testing that would mean looking for local duplication of, or inverted-repeat similarity to, other transcribed loci (proto-MIR origin) before any candidate is promoted.
3. **tRNA CCA FASTA.** It is consistent with the documented rule: 717 sequences have CCA appended and 164 already ended in CCA.

## What this does and does not establish

- **Established:** the GFF3 coordinates are compatible with this genome file. Each feature's own release-assigned sequence occurs exactly at its own coordinates in a single documented orientation. GFF ↔ FASTA ↔ equivalence table ↔ Phasis observation links are internally consistent.
- **Not established:**
  - Biological class: a sequence occurring at an interval does not show it was produced by Phasis-phased, Dicer, MIR, or structural-RNA processing. That rests on the original Phasis window/register, miRScore and Rfam/tRNAscan evidence, which this run does not re-derive.
  - Which locus a multi-placement read came from (14,270 phasiRNA features and most rRNA placements).
  - Whole-genome identity to the originally hashed FASTA.
  - Anything about the Rambo/Phytozome genome (not re-tested here, by design).

## Caveats

- The local genome hash (`52f06632…`) differs from the report's `a40ee358…`. Conclusion accepted 2026-11-05: same B73v5 assembly with different headers/formatting. All 87,606 features and all coordinate-bearing sequences agree with this file, so this is not a blocker.
- Phasis source was inspected in a checkout labelled 2.8.4. The release commit `00fb8c0` is in its history and the parser is unchanged since then.
