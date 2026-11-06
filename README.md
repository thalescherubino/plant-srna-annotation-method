# Plant small-RNA annotation: an evidence-led curation strategy

This repository describes a transferable strategy for building a genome-anchored plant small-RNA annotation. It reconstructs the approach used for maize B73 RefGen_v5 and provides decision points for organisms with sparse reference annotations. It is a **methods and case-study repository**, not an automated end-to-end pipeline.

![Plant small-RNA annotation workflow](figures/workflow.svg)

The vector figure can be downloaded and used in presentations. Its [accessible text description](docs/workflow.md) explains each step and the feedback loops.

## Start here

1. Read the [general workflow and evidence policies](docs/workflow.md).
2. Follow the [software setup guide](docs/setup.md) for miRador, miRScore, Phasis and structural-RNA tools.
3. Use the [resource decision guide](docs/resource_decisions.md) to plan species-specific data collection and identify missing evidence.
4. Consult the [B73v5 case study](case_studies/b73v5/README.md) for the documented tools, versions, policies, and historical scripts.
5. Review the [release checklist](docs/release_checklist.md) before sharing an annotation.

An AI agent or analyst should first inventory the actual assembly, public databases, annotation sources, and available small-RNA libraries for the target species. Evidence availability determines which branches can be completed. Record failed searches and limitations; do not treat a lack of a database hit as proof that a feature is absent.

## Scope and interpretation

Use [miRador](https://github.com/rkweku/miRador) **and** [ShortStack](https://github.com/MikeAxtell/ShortStack) **4.1.0 or newer** to generate plant miRNA candidates, then [miRScore](https://github.com/Aez35/miRScore) to evaluate suitable candidate mature/precursor models using individual small-RNA libraries. Enable ShortStack's de novo miRNA search with `--dn_mirna`; do not use `--nohp`. Use [Phasis](https://github.com/blakemeyerslab/Phasis) **2.8.4 or newer** for *PHAS* and *PHAS*-like locus discovery and classification. Use [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) for tRNAs and [Infernal](https://github.com/EddyRivasLab/infernal) with [Rfam](https://rfam.org/) models, [RNAcentral](https://rnacentral.org/) and genomic support for snRNAs, snoRNAs and rRNAs.

The workflow distinguishes computationally supported and provisional *MIR* models, structural RNA classes, phased loci, and register-supported phasiRNA products. A genomic interval enriched in a particular read length is not by itself a *PHAS* locus, and reads inside a *PHAS* interval are not automatically phasiRNAs. Identical sequences may arise from multiple loci or classes; retain that ambiguity instead of duplicating abundance.

Release the genomic annotation, sequence references, evidence tables, and exclusion decisions together. Genomic sequence and mature-product sequence can differ, such as post-transcriptional tRNA CCA addition; keep their intended uses explicit. Analysis of differential abundance is a separate downstream task.

## Status

The scripts in the case study are historical maize-specific implementations. They are provided as auditable examples, not advertised as species-independent tools. Their input data are not bundled. Consult the original software licenses and source-data terms before redistributing dependencies or reference files.
