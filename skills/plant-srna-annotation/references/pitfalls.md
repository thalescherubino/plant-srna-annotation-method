# Pitfalls and troubleshooting: symptom → cause → fix

[Back to the skill map](../SKILL.md)

Every entry happened on a real run (maize, rice or tomato). Search this page for the symptom you see.

## Library sourcing and processing

| Symptom | Cause | Fix |
|---|---|---|
| Library has reads of 50-150 nt and no 21/24-nt peak after trimming | long-RNA library filed as miRNA-Seq/ncRNA-Seq | exclude; judge text per experiment and require a small-RNA size profile ([02](02_library_sourcing.md)) |
| Same sample appears twice | GEO SuperSeries repeats sub-series samples, or a file was re-deposited | drop SuperSeries copies; checksum files ([02](02_library_sourcing.md)) |
| Counts are decimals | normalized file (RPM/TPM) | exclude; only raw integer counts |
| Only 20-21-nt reads, or only known miRNAs | truncated or miRNA-only table | exclude |
| Converted library has wrong totals | count-first column order read as sequence-first | detect layout per file before converting |
| A false 23-nt peak where 24 nt was expected | over-trimming: an extra 3′ base was cut to "rescue" reads with untemplated 3′ additions | accept a cut only if it at least doubles exact-mapped reads ([03](03_library_processing_qc.md)) |
| Almost no reads map although the size profile looks fine | 4N random ends, inline barcode or a fixed extra base not removed | let the trimming evaluation test 5′/3′ cuts; these give huge mapping jumps |
| 24-nt mapping much lower in some libraries | non-reference cultivar; 24-nt siRNAs come from fast-evolving repeats | keep if they add tissues; record per library; weaker evidence. Ask before swapping libraries |
| Reads > 24 nt map poorly in every library | genome has no organelle contigs; abundant organelle fragments cannot map | expected; not a genotype effect |
| One tool saw different reads from another | tools were fed different files | derive every tool's input from the same frozen tag counts ([03](03_library_processing_qc.md)) |
| miRador used an older library set | its ini pointed at an unfrozen folder | check `libFolder` against the frozen set before running |

## Compute and queue

| Symptom | Cause | Fix |
|---|---|---|
| `nproc` says 1 on a big server | `OMP_NUM_THREADS=1` | use `getconf _NPROCESSORS_ONLN`; check the load |
| Queue did strange things after a script edit | bash reads a running script incrementally | never edit a running script; write a new one |
| `pkill -f pattern` killed your own shell | the pattern matched your command line | kill by PID |
| Step logged `done` but output missing | an earlier command in the step's function failed and the last one succeeded | end each step function with an explicit exit-status check |
| Job died when the laptop slept | run was attached to the session | run in `tmux` or `setsid nohup` with a log |
| Server load far above your share | several heavy jobs at once, or other users | one heavy job at a time; ask before running steps in parallel |
| Untested pandas 3.x in the Phasis env | Phasis does not pin pandas | install with the pip constraints file |

## miRNA prediction and miRScore

| Symptom | Cause | Fix |
|---|---|---|
| ShortStack reports no *MIRNA* loci | `--nohp` was set, or `--dn_mirna` missing | use `--dn_mirna`; never `--nohp` |
| miRador overwrote an earlier result | output folder already existed | always a new output folder |
| miRador fails at `samtools faidx` | samtools not on `PATH` | run inside the environment (`mamba run`) |
| miRScore exits 0 but writes no results | a hairpin has no name-matched mature (miRBase `.1`/`.2` isoforms) | run `hairpinHelper`; count a chunk done only if its results CSV exists |
| hairpinHelper writes no file | no hairpin needed new entries | use the original hairpins ([06](06_mirna_validation_mirscore.md)) |
| miRScore stops: "These miRNAs start at position 1 of the hairpin" | the supplied hairpin starts with the mature, so the 3′ overhang cannot be checked | extend the hairpin from the genome; if it is not found exactly, drop it and log it in `not_scorable.tsv` |
| miRScore treats FASTA libraries as FASTQ | miRScore 0.3.7 reads every `.gz` input as FASTQ | give uncompressed FASTA |
| Precursors > 300 nt fail | plant mode limit | cap hairpins at 300 nt (input builder does this) |
| One miRScore chunk runs for hours with a frozen log | `-rescue` counts reads per library with one `samtools view` call each; very abundant miRNAs (miR166) make each call slow | not a hang: check that `samtools` is using CPU; wait |
| Merged results have missing or duplicated hairpins | chunk results from different splits were mixed | after any input change, re-split and rerun all chunks |
| Many miRador candidates flagged "23/24 nt miRNA" | 24-nt hairpin-derived sRNAs | review separately; real 24-nt miRNAs exist but are rare |

## Phasis

| Symptom | Cause | Fix |
|---|---|---|
| Index rebuilt for every run | each run used its own directory | one run directory per species, runs in series, separate `--outdir` |
| Corrupted cache or odd results | two Phasis runs in one directory at once | never run two at once in one directory |
| 24-nt run takes more than a day | about 100× more clusters pass the pre-filter than at 21 nt; scoring and feature assembly scale with that | expected; schedule other steps meanwhile; report tool bottlenecks to developers |
| Pooled 24-nt run adds thousands of broad loci | dense hc-siRNA clusters produce chance phasing | keep the run; choose a Howell cut-off from the distributions with the project lead |
| phasiRNA sequences on `c` (reverse-strand) rows do not match the library reads | Phasis writes `tag_seq` as aligned, i.e. reference-forward | reverse-complement `c` rows once when you need the sRNA 5′→3′ (Phasis `docs/phasirna_output.md`); never flip `w` rows or flip twice; state the orientation in any release |
| *PHAS* counts disagree between reports | GFF loci vs phased units, or `calls.tsv` rows (locus × library) vs loci | say which unit you count; count distinct `identifier` for loci |

## Structural RNAs and comparisons

| Symptom | Cause | Fix |
|---|---|---|
| No organelle-type rRNA at all | database-first workflow excluded organelle-labelled records | predict genome-wide first |
| Hundreds of hits of one snoRNA family | repeat-derived copies | review copy number; partial hits are exploratory |
| Known RNAcentral ncRNA "missing" | exact full-length matching fails across cultivars (SNPs, indels) | report as no exact match, not as absent; near-identical matching is planned |
| A miRBase entry does not match the genome exactly | the entry comes from another cultivar or assembly | expected for some; log it, do not edit the sequence |
