# Step 10: Build and validate the release

[Back to the skill map](../SKILL.md) · Previous: [09_outputs_and_interpretation.md](09_outputs_and_interpretation.md)

## What a release contains

- Per-class and unified **GFF3** files on the reference genome.
- **FASTA** files per class: genomic sequences, plus mature-product sequences where they differ (spliced tRNAs with CCA, mature miRNAs, phasiRNAs).
- Keep *PHAS*-like phasiRNAs (`<P>_PHAS_like_phasiRNAs.tsv`) in their own files, separate from *PHAS* phasiRNAs, as Phasis's own output does. State the sequence orientation used: Phasis writes `tag_seq` as aligned to the reference; the original read is the reverse complement for `c` (bottom-strand) rows.
- **Evidence tables:** miRScore results, Phasis calls and evidence, structural hit tables, library manifest and QC.
- **Ledgers:** excluded libraries, dropped or unscorable candidates (`not_scorable.tsv`), held-for-review items, each with its reason.
- **Manifests:** software versions and commits, database releases, commands, input and output checksums.
- A **changelog** and a validation report.

## Tiers

Publish canonical, exploratory and organelle-type sets separately, with their decision rules. In this project all tiers serve as truth sets for model validation; none is discarded.

## Validate before sharing

- [ ] Every feature's sequence, re-extracted from the genome by its ID and coordinates, matches its FASTA entry. Check by ID, not just by interval overlap.
- [ ] FASTA alphabet and unique IDs; GFF3 coordinates within sequence lengths; parent-child links resolve.
- [ ] Sequence names match the genome exactly.
- [ ] No class leaks into another file (for example organelle-type rRNA in canonical rRNA).
- [ ] Shared sequences across loci are marked, not double-counted.
- [ ] Checksums recorded; validation re-run after packaging.

## Wording

- "Supported by miRScore on N libraries", not "validated".
- "*PHAS* locus called by Phasis (21 nt)", not "phasiRNA-producing gene with function X".
- State limitations: missing tissues, non-reference genotypes, absent organelle contigs, branches not completed.
- Protect unpublished data and respect the licences of databases and tools.
