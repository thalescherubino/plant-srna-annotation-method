# Step 2: Find and select public sRNA libraries (GEO and SRA)

[Back to the skill map](../SKILL.md) · Previous: [01_reference_genome.md](01_reference_genome.md) · Next: [03_library_processing_qc.md](03_library_processing_qc.md)

## Goal

Up to 100 genuine small-RNA-seq libraries covering as many organs and tissues as possible, each with a recorded reason for being chosen, and every rejected sample listed with its reason.

## Scripts (repository `scripts/libraries/`)

| Order | Script | What it does |
|---|---|---|
| GEO route | `geo_inventory.py` | one row per GEO sample with metadata |
| | `curate_geo_candidates.py` | rule-based labels: cultivar, sRNA library yes/no, genotype class, organ, treatment |
| | `probe_geo_formats.py` | reads the first 64 KB of each processed file and classifies its layout |
| SRA route | `sra_inventory.py` | one row per SRA experiment (SRX), keeping the design text that names the kit |
| both | `select_libraries.py` | applies the exclusions and picks up to N libraries round-robin across organs; writes an exclusion ledger |
| both | `species_profiles.py` | species rules (`--species rice` or `tomato`); add a profile for a new species |

## GEO or SRA?

- **GEO processed files first** if there are enough read-level ones (sequence plus raw count per distinct read). Rice had enough: 100 libraries from 18 series.
- **Raw SRA reads** when processed files are too few or too uniform. Tomato GEO had only 21 read-level samples (14 root, 7 leaf), so all tomato libraries were built from SRA with one uniform pipeline ([03_library_processing_qc.md](03_library_processing_qc.md)). Do not mix half GEO-processed and half SRA-processed libraries when you can avoid it: different upstream trimming makes libraries incomparable.

## Selection policy

- **At most 100 libraries, maximum organ diversity.** Organs are picked round-robin so leaves do not crowd out roots and flowers. Untreated wild type ranks first; then non-deleterious lines (near-isogenic, introgression, hybrids) and treated samples.
- **Genotype does not need to match the reference**, but costs reads under exact matching. Measure it (step 3) and record it per library.

## Pitfalls when sourcing public sRNA data

These all happened. Each now has an explicit rule in `select_libraries.py`.

| Pitfall | What it looks like | Rule |
|---|---|---|
| **Long-RNA libraries under small-RNA labels** | SRA strategy "miRNA-Seq" or "ncRNA-Seq" but the text says lncRNA, circRNA, whole transcriptome, oligo-dT, rRNA-depleted, or reads are ≥ 100 nt. 102 of 1,041 tomato experiments. | Judge text **per experiment** (studies mix library types), then require a small-RNA size profile after trimming. |
| **Wrong protocol** | degradome/PARE, csRNA-seq, TRM-seq, immunoprecipitation (AGO-IP) libraries filed as small RNA | exclude: their read populations are biased by design |
| **Not raw counts** | normalized values (RPM/TPM), count tables of known miRNAs or loci only, reads truncated to 20-21 nt | exclude; only read-level files with integer counts |
| **Duplicates** | GEO SuperSeries repeat their sub-series samples (51 rice libraries); the same file re-deposited under another accession (byte-identical) | drop SuperSeries copies; checksum files and drop identical ones |
| **Column order differs** | some files are COUNT then SEQUENCE | detect layout per file; convert, preserving totals |
| **Pre-filtered files** | uploaders removed rRNA/tRNA or kept only mapped reads | allowed, but flagged; record the filter |
| **Mutants and transgenics** | dcl, rdr, ago mutants, overexpressors, sterile lines | exclude obvious deleterious ones; they change exactly the pathways we annotate |
| **Unclear genotype** | for example "revertant" lines | exclude if the genotype cannot be stated |
| **Pathogen or symbiont-treated tissue** | viroid, virus, oomycete, mycorrhiza | allowed (tomato kept 7); foreign reads simply do not map to the plant genome |
| **Unknown tissue** | metadata missing | allowed but labelled `unknown_tissue` |

## Hand-off to the next step

Show the person running the project the selection summary: counts by organ, genotype group, treatment, study and kit; the exclusion ledger; and any manual overrides. **Get approval before processing**, and ask before swapping libraries later.
