# Workflow and evidence policy

## 1. Define reference and evidence

Choose a versioned genome assembly; record accession, provider, release, checksums, chromosome/contig identifiers, ploidy and organelle scope. Collect available genome annotations, small-RNA libraries and metadata, known ncRNA sequences, and relevant literature. Verify assembly compatibility instead of inferring it from file names. Document library preparation, adapter and UMI handling, quality control, biological context, and replicate structure. Prefer an evidence panel spanning appropriate tissues and developmental stages; report sampling gaps.

## 2. Build candidate evidence in parallel

### miRNA

Collect species-specific mature and precursor annotations where available and track IDs and versions. Run both [miRador](https://github.com/rkweku/miRador) and [ShortStack](https://github.com/MikeAxtell/ShortStack) 4.1.0 or newer on a genome and suitable small-RNA libraries to generate plant miRNA candidates. In ShortStack, enable `--dn_mirna`, optionally supply `--known_miRNAs` when reliable mature sequences are available, and do not set `--nohp`. Preserve independent predictor provenance, reconcile matching candidates and document unique predictions. Project precursors onto the selected assembly, check strand and allele compatibility, and remap product arms. Prepare matching mature and precursor FASTAs and evaluate suitable candidates using [miRScore](https://github.com/Aez35/miRScore) in plant mode with separate, quality-controlled libraries. Inspect pass/fail flags, mature/star support, duplex structure and library-level evidence; carefully review rescue suggestions. Passing computational criteria supports an annotation but does not by itself establish experimental validation or target function. In the maize case miRScore 0.3.7 evaluated reference and published multi-predictor candidates, including repaired genomic models. Merge overlapping same-strand models only under a documented locus rule; keep distinct loci that share a product sequence. Failed or untestable candidates remain exploratory, with reasons. See [setup](setup.md) for environment and input preparation.

### Structural RNA

Use [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) in an appropriate organism mode for genomic tRNAs and review pseudogenes. For snRNAs, snoRNAs and rRNAs, collect species candidates from [RNAcentral](https://rnacentral.org/) where available; map full sequences to the selected genome with [Bowtie 1](https://bowtie-bio.sourceforge.net/index.shtml) under an explicitly chosen mapping policy, and classify families with [Infernal](https://github.com/EddyRivasLab/infernal) `cmscan` and versioned [Rfam](https://rfam.org/) covariance models. The B73 workflow used exact full-length mapping and Rfam gathering thresholds; this strict policy is not universally appropriate when true alleles differ. Record sequence, coordinate, class, source, and ambiguity. A species lacking database entries can use Rfam-guided genome-wide homology searches followed by locus and family review; a related-species match alone is not a validated annotation. Organelle, ambiguous-base, repeat, and pseudogene policies must be explicit. Preserve genomic tRNA and, if justified, a separate mature-read CCA-aware representation.

### Phased small RNAs

Process suitable small-RNA libraries consistently. For new *PHAS* and *PHAS*-like discovery and classification use [Phasis](https://github.com/blakemeyerslab/Phasis) **version 2.8.4 or newer**, running individual-library and, where appropriate, pooled candidate discovery for chosen phase lengths; pooled discovery does not replace biological replicates. Retain *PHAS*-like calls as a distinct, weaker-evidence class, not high-confidence *PHAS* loci. Pin tool version, commit, mapping policy, phase, thresholds and command. Review candidate registers, length, complexity, strand, repeats, locus architecture and library support. Retain phased products only where a tool or validated procedure supplies locus, strand, register and product evidence. Keep locus-only additions locus-only. ShortStack v4 miRNA discovery is a separate branch; its v4 releases do not provide the historical ShortStack 3.8.5 phasing scores.

## 3. Reconcile and review

Compare cross-class genomic overlaps and sequence identities; retain all membership and mapping ambiguity. Define precedence or multi-label policies before release and retain evidence, reviewer decisions, rejected calls, and ties. Maize used specific MIR, structural-RNA and cross-phase precedence rules; these are case-specific decisions, not universal biological truths. Distinguish canonical, exploratory and unsupported classes.

## 4. Build and validate a release

Create per-class and unified genomic GFF3; class-specific genomic and mature sequence FASTAs; product, locus and sequence-equivalence tables; exclusion and uncertainty ledgers; and source and build manifests. Validate sequence alphabet, unique IDs, coordinates and parent-child relationships, assembly/contig compatibility, no unintended class leakage, and checksums. Record software versions, exact commands, settings, inputs, evidence provenance, licenses and review status. Re-run validation after packaging.

## Figure description

The diagram starts with a genome, small-RNA libraries and external resources. Three parallel branches generate miRNA, structural-RNA and *PHAS* evidence. They feed joint cross-class review, which can send ambiguous candidates back for further evidence. The final branch produces a versioned annotation and evidence package. A resource-gap decision applies across the three discovery branches and can yield explicitly provisional outputs.
