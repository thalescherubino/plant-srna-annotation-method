# Decision log

Decisions made while extending the method beyond maize B73v5. Each entry gives the date, the decision and the reason. Later entries supersede earlier ones where they conflict. Historical B73v5 v1 provenance is not rewritten.

## 2026-11-05

| Decision | Reason |
|---|---|
| Validate B73v5 v1 against the local genome feature by feature, using IDs, not only by interval. | An interval match does not prove that a feature carries its own assigned sequence. All 87,606 features passed (see `case_studies/b73v5/validation_2026-11-05.md`). |
| Treat the phasiRNA FASTA and GFF `Name` as **reference-forward**, including minus-strand features; do not edit strands. | Phasis stores the SAM `SEQ` field verbatim (verified in source at commit `00fb8c0`); all 58,123 phasiRNA children match reference-forward. |
| Remove tRNA introns from the mature CCA-aware quantification FASTA; keep them in the GFF3 and the genomic tRNA FASTA. | Mature tRNAs are spliced, so the v1 CCA-aware FASTA could not capture reads crossing the exon–exon junction. Low-score intron calls (13 features) are held for review. |
| Do not allow within-species mismatches to mature miRNAs; revisit cross-species tolerance later. | The 38 exploratory B73 products that differ by 1–2 nt from their locus are source-vs-genome allele artefacts, not paralog variation. |
| Annotate each species on the PlantCAD2 genome FASTAs, even when a newer assembly exists. | The annotations are evaluation truth for PlantCAD2. New genomes are a later phase. |
| Analyse polyploids (e.g. *Coffea arabica*) per subgenome, later. | Homeolog multimapping confounds exact-match evidence. |
| Use miRador as a Python 3.12 port, kept local and not pushed upstream; use conda-forge and bioconda only. | Upstream pins Python 3.7 and Anaconda `defaults`. The port reproduces upstream test outputs byte-for-byte (`docs/setup.md`). |

## 2026-11-06

| Decision | Reason |
|---|---|
| Next species: rice (Nipponbare), then tomato (Heinz 1706), sequentially. | Rice has the best combination of genome and processed public sRNA data. Tomato replaced grapevine; exact genotype matching of libraries is not required. |
| Phasis alone performs *PHAS* analysis, with `-mismat 0`. ShortStack 4 is used for *MIR* prediction and to report how its clusters overlap *PHAS*-like loci; those are reported but not used in later steps. | ShortStack 4 has no phasing analysis. Whether *PHAS*-like and non-*MIR* clusters are untranslated transcripts or degradation products is decided later. |
| Do not use published *PHAS* catalogs as references. | Only loci called and reviewed in the study are annotated. |
| Do not use PmiREN; miRBase (v23.1 at time of use) is the miRNA comparison source. | Consistent with B73v5, which used only miRBase v23, miRScore study S4 predictions and an expanded MIR2275 set (verified in the release). |
| Predict miRNA and structural RNAs ourselves; compare with database records afterwards. | Avoids inheriting database gaps and exclusions. |
| Library policy: ≤ 100 libraries, maximum organ diversity, processed GEO files (pre-filtered allowed and flagged), exclude obvious deleterious mutants and non-standard protocols, record exact file layouts. | See `docs/workflow.md` §1. |
| Keep indica libraries in rice despite lower exact-mapping rates; record rates per library and treat low-mapping libraries as weaker evidence. | Indica libraries provide roots and most vegetative tissues. 24-nt exact mapping: indica 70%, unrecorded cultivar 68%, Nipponbare 94%, other japonica 92%. |
| No gene annotation or organelle genomes are added in this round, and organelle-derived sequence is not flagged in libraries. | Scope decision for the first round. |
| **Structural RNAs for every species: genome-wide prediction (tRNAscan-SE plus Rfam structural-subset `cmscan`), then comparison with RNAcentral.** B73v5 is updated to Version 2 accordingly. | An RNAcentral-first workflow misses organelle-derived rRNA insertions. Rice has 75 bacterial-model rRNA hits on nuclear chromosomes; B73v5 v1 has none, by construction. |
| v2 structural rules: Rfam hits covering ≥ 80% of the model are canonical; partial hits are exploratory; organelle-type rRNA is a separate class outside canonical nuclear rRNA and quantification FASTAs; RNAcentral-only placements stay canonical, flagged; B73v5 v2 absorbs the tRNA spliced patch, with 13 low-score tRNAs held. | See `docs/workflow.md` (structural RNA). |
| Second round (later): sirenRNAs, unknown transcribed regions, degradation-product discrimination. Keep per-locus size, strand, 5′-nucleotide and overlap features now. | Planned research directions, including a possible link between cultivar 24-nt source divergence and maternal-sRNA hybrid-compatibility models. Not acted on yet. |

## 2026-11-07

| Decision | Reason |
|---|---|
| Run every Phasis analysis for a species (21/24 nt, individual/pooled) in **one run directory**, serially, with a distinct `--outdir` per analysis. | Phasis reuses the genome index, processed libraries and its `phasis.mem` cache within a run directory; the index is large and slow to build. Running several analyses in one directory at once is unsafe, but Phasis is expensive enough that it should not be run in parallel anyway. Supersedes the earlier one-directory-per-run setup used for rice. |
| miRScore 0.3.7 in plant mode validates every miRNA set on the same libraries: ShortStack 4 de novo, miRador, and miRBase entries. miRBase is a tested set, not a positive control. `-rescue` is always on. | miRBase may contain false positives. Rescue recovers true MIRNAs whose annotated duplex is not the dominant one. |
| miRScore hairpins for predictor sets: miR/star duplex span ± 20 nt, capped so no hairpin exceeds 300 nt. miRBase hairpins are kept as deposited; those whose mature sits at a hairpin end are extended under the same cap. miRBase sets go through `hairpinHelper` first. | miRScore's plant mode fails precursors over 300 nt, and precision depends on the hairpin span, so both predictors get the same rule. miRScore exits 0 without results when a hairpin has no name-matched mature (miRBase `.1`/`.2` isoforms). |
| miRScore libraries are expanded tag counts in **uncompressed** FASTA. Large sets may be split into hairpin chunks run in parallel and then merged. | miRScore 0.3.7 maps every `.gz` input as FASTQ. Each hairpin is scored independently and bowtie `-a` reports every alignment, so chunking does not change results. |
| The current method is "Version 2" of the workflow. Version 1 is the maize B73v5 build from the Meyers lab; maize will be redone with Version 2. | Naming for presentations and the record. |
| Run `phasis-compare` on the pooled and individual runs for each phase; the combined output is the final *PHAS* set (default `--prefer auto`, which keeps the non-pooled record for shared loci). | Phasis docs, "Combine Completed Runs". Rice 21-*PHAS*: 273 pooled + 1,000 individual → 1,011 loci (808 shared, 192 individual-only, 11 pooled-only). |
| Rename the organelle class to **organelle-type ncRNAs** (bacterial-like rRNAs and tRNAs). Whether a genome includes organelle contigs varies, so it is never assumed. | Organelle-type sequences can be genuine organelle genes (when organelle contigs are in the assembly) or nuclear copies of organelle DNA (NUPT/NUMT). Organelles also have their own tRNAs, so the class is not limited to rRNA. |
| All annotation sets (canonical, exploratory, organelle-type) are truth sets for validating PlantCAD2; none is treated as lesser. "Quantification set" is no longer the goal. | The project goal changed from quantification to model validation. |

## Planned improvements (not yet implemented)

| Topic | Plan |
|---|---|
| RNAcentral comparison | The current exact, full-length matching is too strict across varieties: a cultivar may carry SNPs or indels in a conserved ncRNA. Replace it with near-identical alignment (identity and coverage thresholds) that reports the differences. Variant analysis comes later. |
| Organelle-type tRNAs | When organelle contigs are present, run tRNAscan-SE in organellar mode on them. For nuclear tRNAs, flag organelle-type (bacterial-like) candidates by best match to organelle tRNAs. |
| hc-siRNA loci | ShortStack 4 gives a DicerCall, not an hc-siRNA label. Treat 24-nt-dominant clusters as candidate heterochromatic siRNA loci, as studies usually do (e.g. Lunardon et al. 2020, Genome Res, 47 plants grouped by dominant size). Support each call with overlap with transposons and repeats, a 5′-A bias, and DNA methylation at the locus (PacBio methylation data). |
| Degradation fragments vs bona fide loci | First candidates are ShortStack DicerCall "N" clusters (broad, non-Dicer read-size spread). Add read-size distribution shape, strand bias and 5′-nucleotide profiles as evidence. |
| Compute planning | About 2.5 days per species on bl02 at 32 cores (rice, 374 Mb). All 11 species take about 25–55 days, depending on genome size (Phasis 24-nt runs dominate) and on how much QC-passed GEO sRNA-seq data exists. |
