# Workflow and evidence policy

## 1. Define reference and evidence

Choose a versioned genome assembly; record accession, provider, release, checksums, chromosome/contig identifiers, ploidy and organelle scope. Collect available genome annotations, small-RNA libraries and metadata, known ncRNA sequences, and relevant literature. Verify assembly compatibility instead of inferring it from file names. Document library preparation, adapter and UMI handling, quality control, biological context, and replicate structure. Prefer an evidence panel spanning appropriate tissues and developmental stages; report sampling gaps. Polyploids should be analysed by splitting subgenomes and treating each as a separate entity, because homeolog multimapping otherwise confounds exact-match evidence.

### Public small-RNA libraries

The scripts in [`scripts/libraries/`](../scripts/libraries/) implement the default library policy, adopted 2026-11-06:

- **Cap and diversity:** at most 100 libraries per species, chosen to maximize organ and tissue diversity.
- **Processed data only:** use processed GEO files (pre-filtered files are allowed but flagged). Never use count tables of miRNAs, genes or loci.
- **Exclusions:**
  - obvious deleterious mutants, transgenics and sterile lines (non-deleterious near-isogenic, introgression and hybrid lines are allowed, ranked after wild type);
  - non-standard protocols (degradome/PARE, csRNA-seq, TRM-seq, immunoprecipitation);
  - normalized values;
  - libraries truncated to miRNA lengths;
  - duplicate deposits (SuperSeries repeats and byte-identical files under different accessions).
- **Format conversion:** record each file's exact layout before use. Different aligners treat inputs differently, so every tool receives the same reads, derived from one standardized tag-count file per library with abundances preserved.
- **Genotype mismatch:** library genotypes need not match the reference. Measure the cost with the exact-mapping check ([`scripts/qc/`](../scripts/qc/)) and record it per library. In rice, about 30% of indica 24-nt reads did not map exactly to Nipponbare, against 6–8% for japonica.

## 2. Build candidate evidence in parallel

### miRNA

Collect species-specific mature and precursor annotations where available and track IDs and versions. Use them to compare with our own predictions, not as the starting point. miRBase is the comparison source; PmiREN was not used for B73v5 and is not used by default. Run both [miRador](https://github.com/rkweku/miRador) and [ShortStack](https://github.com/MikeAxtell/ShortStack) 4.1.0 or newer on a genome and suitable small-RNA libraries to generate plant miRNA candidates. In ShortStack, enable `--dn_mirna`, optionally supply `--known_miRNAs` when reliable mature sequences are available, and do not set `--nohp`. Preserve independent predictor provenance, reconcile matching candidates and document unique predictions. Project precursors onto the selected assembly, check strand and allele compatibility, and remap product arms. Prepare matching mature and precursor FASTAs and evaluate suitable candidates using [miRScore](https://github.com/Aez35/miRScore) in plant mode with separate, quality-controlled libraries. Inspect pass/fail flags, mature/star support, duplex structure and library-level evidence; carefully review rescue suggestions. Passing computational criteria supports an annotation but does not by itself establish experimental validation or target function. In the maize case miRScore 0.3.7 evaluated reference and published multi-predictor candidates, including repaired genomic models. Merge overlapping same-strand models only under a documented locus rule; keep distinct loci that share a product sequence. Failed or untestable candidates remain exploratory, with reasons. See [setup](setup.md) for environment and input preparation.

### Structural RNA

**Policy, adopted 2026-11-06 for every species, including the B73v5 Version 2 update: predict genome-wide first, then compare with RNAcentral.** An RNAcentral-first workflow cannot find features missing from, or excluded from, the database. The B73v5 v1 structural branch started from RNAcentral and excluded organelle-labelled records, so it contains no organelle-derived rRNA insertions. A genome-wide Rfam scan of rice found 75 bacterial-model rRNA hits, all on assembled nuclear chromosomes, consistent with chloroplast or mitochondrial insertions (NUPTs/NUMTs).

1. **tRNAs:** [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) in eukaryotic mode on the whole genome. Keep scores, isotypes and pseudogene calls, and flag (not remove) predictions failing EukHighConfidenceFilter.
2. **rRNA, snRNA and snoRNA:** genome-wide [Infernal](https://github.com/EddyRivasLab/infernal) `cmscan --cut_ga --rfam --nohmmonly --fmt 2 --clanin`, with Z = 2 × genome Mb.
   - Models: the versioned [Rfam](https://rfam.org/) structural subset (rRNA, spliceosomal snRNA, C/D and H/ACA snoRNA, scaRNA; 827 families in Rfam 15.1).
   - Keep non-overlapping hits after clan competition.
   - See [`scripts/structural/`](../scripts/structural/).
3. **Comparison with [RNAcentral](https://rnacentral.org/):** compare afterwards, using typed records mapped full-length with [Bowtie 1](https://bowtie-bio.sourceforge.net/index.shtml) under an explicit mapping policy (exact for B73), plus the release's genome-coordinate file where available.
4. **Rules:**

| Case | Treatment |
|---|---|
| Rfam hit covering ≥ 80% of its model, nuclear-type best model | canonical; RNAcentral support and model coverage as attributes |
| Rfam hit covering < 80% of its model | exploratory (`partial_rfam_hit`); fragments are often repeat-derived |
| Best model bacterial/archaeal LSU/SSU rRNA | separate class `organelle_type_rRNA` (putative NUPT/NUMT); not canonical nuclear rRNA; not in quantification FASTAs |
| Existing RNAcentral-based placement without an Rfam hit | kept canonical, flagged `rfam_genome_scan=no_hit` for review |

Record sequence, coordinate, class, source, and ambiguity. A related-species match alone is not a validated annotation. Ambiguous-base, repeat and pseudogene policies must be explicit. For tRNAs, keep the genomic (coordinate) representation and a separate mature-read quantification FASTA:
- one terminal CCA is added only when not already encoded;
- introns are removed from canonically spliced tRNAs;
- identical mature sequences are collapsed;
- low-score intron calls are held for review.

### Phased small RNAs

Process suitable small-RNA libraries consistently. For new *PHAS* and *PHAS*-like discovery and classification use [Phasis](https://github.com/blakemeyerslab/Phasis) **version 2.8.4 or newer**, running individual-library and, where appropriate, pooled candidate discovery for chosen phase lengths; pooled discovery does not replace biological replicates. Retain *PHAS*-like calls as a distinct, weaker-evidence class, not high-confidence *PHAS* loci. Pin tool version, commit, mapping policy, phase, thresholds and command. Review candidate registers, length, complexity, strand, repeats, locus architecture and library support. Retain phased products only where a tool or validated procedure supplies locus, strand, register and product evidence. Keep locus-only additions locus-only. ShortStack v4 releases do not provide the historical ShortStack 3.8.5 phasing scores, so **Phasis alone performs the *PHAS* analysis.** Run it with exact mapping (`-mismat 0`, the default) and record `-maxhits`. **Run all Phasis analyses for a species (21 and 24 nt, individual and pooled) in one run directory**, one after another, separating results with `--outdir` (for example `21_results_individual`, `24_results_pooled`). Phasis is designed to reuse its genome index, processed libraries and `phasis.mem` cache within a run directory, and building the index is the most expensive step. Simultaneous runs in one directory are unsafe, so use a second run directory only when two analyses must run at the same time, which repeats that cost. ShortStack 4 is used for *MIR* prediction and to report how its other small-RNA clusters overlap Phasis *PHAS*-like loci. Those *PHAS*-like loci and non-*MIR* ShortStack clusters are reported but not used in later steps until reviewed. Published *PHAS* catalogs are not used as references; only loci called in the study are annotated.

## 3. Reconcile and review

Compare cross-class genomic overlaps and sequence identities; retain all membership and mapping ambiguity. Define precedence or multi-label policies before release and retain evidence, reviewer decisions, rejected calls, and ties. Maize used specific MIR, structural-RNA and cross-phase precedence rules; these are case-specific decisions, not universal biological truths. Distinguish canonical, exploratory and unsupported classes.

## 4. Build and validate a release

Create per-class and unified genomic GFF3; class-specific genomic and mature sequence FASTAs; product, locus and sequence-equivalence tables; exclusion and uncertainty ledgers; and source and build manifests. Validate sequence alphabet, unique IDs, coordinates and parent-child relationships, assembly/contig compatibility, no unintended class leakage, and checksums. Record software versions, exact commands, settings, inputs, evidence provenance, licenses and review status. Re-run validation after packaging.

## Figure description

The diagram starts with a genome, small-RNA libraries and external resources. Three parallel branches generate miRNA, structural-RNA and *PHAS* evidence. They feed joint cross-class review, which can send ambiguous candidates back for further evidence. The final branch produces a versioned annotation and evidence package. A resource-gap decision applies across the three discovery branches and can yield explicitly provisional outputs.
