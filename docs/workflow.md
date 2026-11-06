# Workflow and evidence policy

## 1. Define reference and evidence

Choose a versioned genome assembly; record accession, provider, release, checksums, chromosome/contig identifiers, ploidy and organelle scope. Collect available genome annotations, small-RNA libraries and metadata, known ncRNA sequences, and relevant literature. Verify assembly compatibility instead of inferring it from file names. Document library preparation, adapter and UMI handling, quality control, biological context, and replicate structure. Prefer an evidence panel spanning appropriate tissues and developmental stages; report sampling gaps.

## 2. Build candidate evidence in parallel

### miRNA

Collect species-specific mature and precursor annotations where available, track IDs and versions, and generate de novo candidates when appropriate. Project candidate precursors onto the selected assembly, check strand and allele compatibility, remap product arms, and evaluate hairpin and library-derived evidence using a suitable plant miRNA validation framework. In the maize case miRScore 0.3.7 evaluated reference and de novo candidates, including repaired genomic models. Merge overlapping same-strand models only under a documented locus rule; keep distinct loci that share a product sequence. Failed or untestable candidates remain exploratory, with reasons.

### Structural RNA

Predict tRNAs in the genome and review pseudogenes. Compile available rRNA, snRNA and snoRNA candidates; verify genomic support and family evidence using appropriate sequence and covariance-model methods. Record sequence, coordinate, class, source, and ambiguity. A species lacking database entries can use genome-wide prediction or homology-based candidates followed by family-model and locus review; a related-species match alone is not a validated annotation. Organelle, ambiguous-base, repeat, and pseudogene policies must be explicit. Preserve genomic tRNA and, if justified, a separate mature-read CCA-aware representation.

### Phased small RNAs

Process suitable small-RNA libraries consistently. Discover phased candidates separately by library and by pooling for chosen phase lengths; pooled discovery does not replace biological replicates. Pin tool version, mapping policy, phase, thresholds and command. Compare complementary callers only after auditing assembly, inputs, multimapping, boundaries, scores and criteria. Review candidate registers, length, complexity, strand, repeats, locus architecture and library support. Do not infer false negatives simply from unequal call counts. Retain phased products only where a tool or validated procedure supplies locus, strand, register and product evidence. Keep locus-only additions locus-only.

## 3. Reconcile and review

Compare cross-class genomic overlaps and sequence identities; retain all membership and mapping ambiguity. Define precedence or multi-label policies before release and retain evidence, reviewer decisions, rejected calls, and ties. Maize used specific MIR, structural-RNA and cross-phase precedence rules; these are case-specific decisions, not universal biological truths. Distinguish canonical, exploratory and unsupported classes.

## 4. Build and validate a release

Create per-class and unified genomic GFF3; class-specific genomic and mature sequence FASTAs; product, locus and sequence-equivalence tables; exclusion and uncertainty ledgers; and source and build manifests. Validate sequence alphabet, unique IDs, coordinates and parent-child relationships, assembly/contig compatibility, no unintended class leakage, and checksums. Record software versions, exact commands, settings, inputs, evidence provenance, licenses and review status. Re-run validation after packaging.

## Figure description

The diagram starts with a genome, small-RNA libraries and external resources. Three parallel branches generate miRNA, structural-RNA and PHAS evidence. They feed joint cross-class review, which can send ambiguous candidates back for further evidence. The final branch produces a versioned annotation and evidence package. A resource-gap decision applies across the three discovery branches and can yield explicitly provisional outputs.
