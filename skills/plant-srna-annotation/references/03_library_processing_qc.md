# Step 3: Trim, collapse and QC the libraries

[Back to the skill map](../SKILL.md) · Previous: [02_library_sourcing.md](02_library_sourcing.md) · Next: [04_compute_and_queue.md](04_compute_and_queue.md)

## Goal

One tag-count file per library (`SEQUENCE<TAB>COUNT`, DNA alphabet, no header, 18-34 nt) plus a QC report, then a **frozen** library set that every tool reads.

## Scripts (repository `scripts/libraries/` and `scripts/qc/`)

| Script | Use |
|---|---|
| `inspect_and_standardize_libraries.py` | GEO route: parse any layout, report depth, length profile, 5′ base, adapter traces; write a tag-count file only when counts are raw integers |
| `sra_to_tagcount.py` | SRA route: download read 1, choose trimming from evidence, trim, keep 18-34 nt, collapse; JSON report with QC flags |
| `library_qc_summary.py` | one QC table for the whole selection, with review-only soft flags |
| `freeze_library_set.py` | freeze the approved set: tag-count links, manifest with provenance, ShortStack condensed FASTA, miRador library links |
| `qc/exact_mapping_qc.py`, `qc/plot_mapping_by_length.py` | read-weighted fraction of reads mapping exactly, overall and per length, per library and per cultivar group |

## How trimming is chosen (SRA route)

Kits differ: 3′ adapter sequence, an extra fixed 5′ or 3′ base, random 4N ends, inline barcodes, UMIs, poly-A. Instead of trusting metadata, `sra_to_tagcount.py` **tests** options on a 100,000-read sample:

1. Detect known adapters by k-mer counts (TruSeq small RNA, NEBNext/universal, QIAseq, BGI/MGI, others), or a de novo over-represented 12-mer.
2. Try every combination of adapter × extra 5′ cut (0-4 nt) × extra 3′ cut (0-4 nt). For each, count 18-34-nt reads that map **exactly** to the genome.
3. Start from "no extra cut". Accept a further cut **only if it at least doubles** the mapped reads (JUMP = 2.0), and take the smallest cut within 10% of the best (NEAR = 0.90).

**Why the doubling rule:** a real artefact (4N ends, a fixed extra base) makes most reads unmappable, so removing it gives a huge jump (about 3-4× for a fixed base, more than 1,000× for 4N ends). Cutting one more base *always* gains 20-40% because it removes the cell's own untemplated 3′ additions. Following those small gains turned the 24-nt peak into a false 23-nt peak in 17 tomato libraries under the first rule. The gains are bimodal, so the threshold sits in the gap.

## QC: what to look at per library

| Metric | Healthy | Worry when |
|---|---|---|
| Size profile | peaks at 21 and/or 24 nt | no 21/22 or 23/24 peak; broad smear; peak at 23 nt (over-trimming) |
| 5′U among 21-nt reads | ≥ 0.40 | lower: miRNAs missing or contamination |
| 5′A among 24-nt reads | ≥ 0.40 | lower: few hc-siRNAs, or another class dominates |
| Exact mapping (18-30 nt) | tomato median 88.5%; rice Nipponbare 94% at 24 nt | under 30% of 20-24-nt reads is a hard fail; 40-60% is weaker evidence |
| Top sequence share | a few % | ≥ 20%: one rRNA/tRNA fragment or one miRNA dominates and skews every profile |
| Reads kept | millions | very few reads after trimming |

**Exact mapping depends on genotype.** In rice, 24-nt reads of indica cultivars mapped exactly to the japonica reference at about 70% against 94% for Nipponbare; 21-nt reads lose much less. 24-nt siRNAs come from fast-evolving repeats, so cultivar differences hit them hardest. Keep such libraries if they add tissues, record the rate, and treat them as weaker evidence. Reads longer than 24 nt map poorly when the genome lacks organelle contigs; that is not a genotype effect.

## Freeze

Once approved, freeze the set (`freeze_library_set.py`). From then on every tool reads only the frozen folder:

- `tag/<LIB>.tag`: Phasis input;
- `shortstack_condensed/<LIB>.fa.gz`: ShortStack 4.1 condensed FASTA (`>LIB_Cd<n>_<count>`);
- `mirador_libs/`: miRador inputs;
- later, `expanded_fasta/<LIB>.fa`: one FASTA record per read, **uncompressed**, for miRScore ([06_mirna_validation_mirscore.md](06_mirna_validation_mirscore.md)).

A rice run accidentally pointed miRador at an earlier, unfrozen selection. The queue now refuses to start miRador unless its configuration points at the frozen set.
