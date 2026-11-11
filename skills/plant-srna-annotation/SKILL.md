---
name: plant-srna-annotation
description: Plan, run, troubleshoot and interpret a genome-wide plant small-RNA (sRNA) annotation with the "Version 2" evidence-led workflow - public sRNA library sourcing from GEO/SRA, trimming and QC, miRNA prediction (ShortStack 4, miRador), miRNA validation (miRScore), phased siRNA/PHAS loci (Phasis), structural RNAs (tRNAscan-SE, Infernal/Rfam), and release. Use this skill whenever someone works on plant sRNA-seq, miRNAs, MIR genes, phasiRNAs, PHAS loci, hc-siRNAs, sRNA library selection or trimming, miRScore/Phasis/ShortStack/miRador runs or their outputs, or asks what such results mean, even if they do not name the workflow.
---

# Plant small-RNA annotation (Version 2 workflow)

This skill lets someone with no prior sRNA experience run the annotation, understand each step and avoid the mistakes that cost us days on rice and tomato. Read this page first; it is a map. Open a sub-file only when you reach that step, so you load only what you need.

## What the workflow produces

For one plant genome, a set of genomic annotations backed by evidence from about 100 public small-RNA libraries:

| Class | Evidence | Main tool |
|---|---|---|
| miRNA genes (*MIR*) and their mature miRNAs | predicted de novo, then scored on every library | ShortStack 4, miRador, miRScore |
| *PHAS* loci and their phasiRNAs (21-nt and 24-nt) | phased read patterns | Phasis |
| Structural ncRNAs: tRNA, rRNA, snRNA, snoRNA, organelle-type ncRNA | genome-wide covariance models | tRNAscan-SE, Infernal + Rfam |

All three branches start from the same genome and the same tag-count libraries. Database records (miRBase, RNAcentral) are compared **after** prediction, never used as the starting point, because databases have gaps and their own errors.

New to the vocabulary? Read [glossary.md](references/glossary.md) first (10 minutes).

## The workflow at a glance

```
genome  ──► 1. reference prep ──────────────────────────────┐
public sRNA-seq ─► 2. find & select ─► 3. trim, collapse, QC ─┤
                                                            ▼
             ┌──────────── same tag-count libraries ─────────────┐
             ▼                       ▼                           ▼
   5. miRNA prediction     7. PHAS loci (Phasis)       8. structural RNAs
      (ShortStack, miRador)   21/24 nt, individual        (tRNAscan-SE,
             ▼                and pooled, then compare     Rfam cmscan)
   6. miRNA scoring (miRScore)       ▼                           ▼
             └──────────────► 9. read & review results ◄─────────┘
                                     ▼
                              10. release package
```

## Step guide: which file to read

| Step | Read | Typical time (tomato, 800 Mb genome, 100 libraries, 32 cores) |
|---|---|---|
| 1. Reference genome | [01_reference_genome.md](references/01_reference_genome.md) | minutes |
| 2. Find and select public libraries | [02_library_sourcing.md](references/02_library_sourcing.md) | 1-2 days, mostly review |
| 3. Trim, collapse, QC | [03_library_processing_qc.md](references/03_library_processing_qc.md) | hours (download-bound) |
| 4. Environments, compute, the run queue | [04_compute_and_queue.md](references/04_compute_and_queue.md) | setup 1-2 h |
| 5. miRNA prediction | [05_mirna_prediction.md](references/05_mirna_prediction.md) | ShortStack 2.6-3 h per run; miRador about 1 h |
| 6. miRNA validation with miRScore | [06_mirna_validation_mirscore.md](references/06_mirna_validation_mirscore.md) | 3-5 h per set |
| 7. *PHAS* loci with Phasis | [07_phas_phasis.md](references/07_phas_phasis.md) | 21 nt about 4 h; 24 nt more than a day |
| 8. Structural RNAs | [08_structural_rna.md](references/08_structural_rna.md) | about 7 h |
| 9. Reading the outputs | [09_outputs_and_interpretation.md](references/09_outputs_and_interpretation.md) | |
| 10. Release | [10_release.md](references/10_release.md) | |
| Something went wrong | [pitfalls.md](references/pitfalls.md): symptom → cause → fix | |

## Rules that hold at every step

These exist because breaking them silently corrupts results. Each is explained in the step file where it applies.

1. **One read set for every tool.** Each library becomes one tag-count file (`SEQUENCE<TAB>COUNT`, DNA letters, no header). Every tool's input is derived from that file, so all branches see identical reads.
2. **Exact matching.** Reads must match the genome with zero mismatches (Phasis `-mismat 0`, miRScore 0 mismatches). Non-reference cultivars therefore lose reads; measure the loss, do not hide it.
3. **Only small-RNA-seq.** Many SRA "miRNA-Seq" entries are long-RNA libraries. Check every library's size profile after trimming.
4. **Predict first, compare later.** miRBase and RNAcentral are comparison sets. miRBase is scored like any other set, not trusted as a positive control.
5. **Never modify source files.** Work on copies of the genome and downloads; record checksums of both.
6. **Record every decision with its reason and date.** Exclusions, dropped candidates and failed steps go in a log, never just deleted. Someone will ask why a library or a miRNA is missing.
7. **A computational pass is not experimental validation.** Say "supported by miRScore", not "validated".

## How an agent should work with this skill

- Before any step, inventory what actually exists (genome version, available libraries, database entries). Missing evidence changes which branches can be completed; report the gaps.
- Check the server load before choosing thread counts, and run long jobs detached with a log and a status file ([04_compute_and_queue.md](references/04_compute_and_queue.md)).
- When a step fails, read its log first, then [pitfalls.md](references/pitfalls.md). Fix the pipeline (inputs, wrapper scripts, parameters); do not patch the third-party tools themselves, and write down any tool problem for its developers instead.
- Ask the person running the project before: swapping or dropping libraries, deleting large outputs, choosing score cut-offs, or running heavy jobs in parallel.
- When reporting, say plainly what ran, what failed and what was skipped. Label numbers "preliminary" until reviewed.

## Explaining results to non-specialists

When someone asks what a result means, answer in this order: what the file contains, what one row is, which columns matter, what a typical value looks like, and what it does **not** prove. [09_outputs_and_interpretation.md](references/09_outputs_and_interpretation.md) has this for every output, with rice and tomato numbers for comparison.

## Project defaults you may need to change

These are the settings of the project this skill was written for. Another group can change them, but should record the change:

| Setting | Default here | Why |
|---|---|---|
| Libraries per species | at most 100, maximum organ diversity | breadth of evidence for a reasonable compute cost |
| CPU budget | 32 cores per heavy job, one heavy job at a time (64 for selected Phasis runs on an idle server) | shared server |
| Genome | the PlantCAD2 genome FASTAs | annotations are truth sets for that model |
| Tool versions | Phasis 2.8.4 (commit 90a3a0c), ShortStack 4.1.2, miRScore 0.3.7, miRador Python 3.12 port, tRNAscan-SE 2.0.13, Infernal 1.1.4, Bowtie 1.3.1, Rfam 15.1, miRBase 23.1 | reproducibility |
| Polyploids | split subgenomes and treat each as a genome | homeolog multimapping breaks exact-match evidence |

The method's dated decision log lives in the repository at `docs/decisions.md`; the reusable scripts referenced in the step files are in `scripts/`.
