# Reusable scripts

These scripts are reusable across species. They are separate from the historical, maize-specific builders in `case_studies/b73v5/scripts/`. They were written and tested on Linux while annotating rice (2026-11-06) and take paths as arguments. Raw data and large outputs are never stored in this repository.

Python scripts use only the standard library, except the plot (matplotlib) and the steps that call external tools (Bowtie 1, Infernal, tRNAscan-SE). Each script's docstring states its inputs and outputs.

## 1. Public small-RNA libraries (`libraries/`)

Run in this order:

| Script | Purpose |
|---|---|
| `geo_inventory.py` | Query GEO DataSets for a species: series metadata plus one row per sample from the series matrices. `--term` gives a custom query; `--reparse` rebuilds the sample table from downloaded matrices. |
| `curate_geo_candidates.py` | Rule-based annotation of each sample, using explicit regular expressions so decisions can be reviewed: subspecies/cultivar, small-RNA library yes/no, genotype class (wild type, mutant/transgenic, special line), normalized organ (explicit tissue field first), treatment class. Cultivar patterns are rice-specific and must be adapted per species. |
| `probe_geo_formats.py` | Stream the first 64 KB of each supplementary file and classify its layout: SEQ<TAB>COUNT, COUNT<TAB>SEQ, collapsed FASTA, reads, FASTQ, or not read-level. |
| `select_libraries.py` | Select up to N libraries with maximum organ diversity (round-robin across organs). It excludes non-read-level files, mutants/transgenics, deleterious lines, non-standard protocols (degradome/PARE, csRNA, TRM, IP), SuperSeries duplicates, byte-identical deposits and normalized/truncated libraries, and writes an exclusion ledger with a reason per sample. |
| `download_geo_libraries.py` | Download the chosen file per library; resumable; checksummed manifest. |
| `inspect_and_standardize_libraries.py` | Parse every layout, report sanity statistics (depth, integer counts, length distribution, 5′ nucleotide, adapter traces) and write a standard Phasis tag-count file (SEQ<TAB>COUNT, DNA, no header) only when counts are raw integers. |
| `freeze_library_set.py` | Freeze a selected set: links, manifest with provenance, ShortStack 4.1 condensed FASTA (`>LIB_Cd<n>_<count>`) and miRador library links derived from the same tag counts. |

## 2. Structural RNAs (`structural/`)

The policy for every species is to predict genome-wide first, then compare with RNAcentral. See `docs/workflow.md`.

| Script | Purpose |
|---|---|
| `make_rfam_structural_subset.py` | Build the Rfam structural subset (rRNA, snRNA, snoRNA, scaRNA; 827 families in Rfam 15.1) and a matching clan file. |
| `run_structural_scan.sh` | tRNAscan-SE `-E` with EukHighConfidenceFilter, then genome-wide `cmscan --cut_ga --rfam --nohmmonly --clanin`, with Z = 2 × genome Mb. |

## 3. Quality control (`qc/`)

| Script | Purpose |
|---|---|
| `exact_mapping_qc.py` | Per library, the read-weighted fraction of 18–30-nt reads mapping with zero mismatches (Bowtie 1 `-v 0 -k 1`), overall and per length. This quantifies reads lost to exact-match policies, for example for non-reference cultivars. |
| `plot_mapping_by_length.py` | Exact-mapping rate by read length per cultivar group: median line with interquartile band. Group labels default to the rice groups. |
