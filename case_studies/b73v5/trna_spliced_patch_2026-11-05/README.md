# Patch: spliced mature tRNA quantification FASTA for B73v5 curated sRNA reference v1

Created 2026-11-05. This is an add-on to the frozen release `B73v5_curated_sRNA_reference_v1_2026-08-27/`, which is **not modified**.
Built by `scripts/build_tRNA_spliced_patch.py` (new code, not original release code); the log is in `build.log`.

## Problem

In v1, the 50 multi-exon tRNA features (39 genomic reference sequences) appear in `fasta/tRNA.genomic.fa` as their **full genomic span, intron included**.
`fasta/tRNA.CCA_aware_quantification.fa` was derived from it, so it carries the intron too.
That is correct for a genomic/coordinate representation. It is **wrong for a mature-tRNA quantification reference**: introns are removed during pre-tRNA splicing, so reads from mature tRNAs that span the exon–exon junction cannot map exactly.

## What is kept and what changes

| File | Status |
|---|---|
| v1 `annotation/*.gff3` (tRNA + exon children) | **unchanged**; exon structure already encodes the intron |
| v1 `fasta/tRNA.genomic.fa` | **unchanged**; remains the coordinate-consistent genomic (pre-tRNA) sequence |
| `annotation/B73v5_sRNA_reference_v1.tRNA_introns.sidecar.gff3` | **new**: 50 explicit `intron` features (Parent = v1 tRNA ID), derived from gaps between v1 exons |
| `fasta/tRNA.CCA_aware_quantification.spliced.fa` | **new**: replaces the v1 CCA-aware FASTA **for quantification** |
| `tables/tRNA_mature_sequence_equivalence.tsv` | **new**: quantification record → genomic reference IDs → GFF tRNA feature IDs |
| `tables/tRNA_intron_patch.tsv` | **new**: per multi-exon tRNA feature, with introns, scores, genomic and spliced sequences, and the decision |

## Rules applied

1. Mature sequence = v1 GFF3 exons concatenated in transcript orientation, taken from the B73v5 genome.
2. One terminal CCA is appended only if the spliced sequence does not already end in CCA (the v1 rule).
3. Identical mature sequences are collapsed into one record. The 35 spliced genomic references differ only in their introns, so they collapse to **5** mature sequences, for example 13 Tyr-GTA references → `…CCA_SPL_0001`. No spliced sequence equals an existing intron-less tRNA.
4. A tRNA is spliced only if:
   - its intron is canonical, meaning the anticodon ends exactly 1 nt before the intron (positions 37/38; true for all 50), and
   - its tRNAscan-SE score is ≥ 50.

   Others keep their v1 record, labelled `held_for_review_v1_unchanged`.

## Result

| | Genomic refs | GFF features |
|---|---:|---:|
| Multi-exon tRNAs | 39 | 50 |
| Spliced (Tyr-GTA ×13 refs, Met-CAT ×21, Val-CAC ×1) | 35 | 37 |
| Held for review | 4 | 13 |
| Records in quantification FASTA: v1 → patch | 881 → **851** | |

## Held for review (decision needed)

- **Leu-TAA on unplaced scaffolds:** `zma_B73v5_tRNA_00803` (10 copies) and `…00804` (1 copy).
  - tRNAscan score 31–36, against about 70 for the five intron-less chromosomal Leu-TAA genes (83 nt).
  - Splicing would give a 70–71-nt "Leu" tRNA. That lacks the long variable arm of type-II Leu tRNAs, so the 50-nt "intron" call is doubtful.
  - These may be tRNA-derived repeats or organelle-derived insertions. The assembly has no organelle contigs, so this could not be tested here.
- **Val-TAC** `9.trna62` (score 37.3) and **Cys-ACA** `6.trna74` (score 30.6): canonical intron position, but low scores and isotypes not usually intron-containing in plants.

Options are to splice them (`--splice-all`), keep them as in v1 (current), or move them to an exploratory/low-confidence tRNA tier.

## Use

- **Quantification:** use `tRNA.CCA_aware_quantification.spliced.fa`. Aggregate counts to tRNA family or isodecoder through `tRNA_mature_sequence_equivalence.tsv`.
- **Coordinates:** keep using the v1 GFF3 and `tRNA.genomic.fa`.
- **Pre-tRNA reads:** reads containing intron sequence (pre-tRNA, 5′ leader/intron fragments) will no longer map to spliced records. Add `tRNA.genomic.fa` as a separate pre-tRNA category if those are of interest.
