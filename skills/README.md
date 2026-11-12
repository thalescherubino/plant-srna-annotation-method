# Skills

AI-agent skills that package this method for bioinformaticians, biologists and AI agents: how to run it, how to read the results and how to avoid pitfalls.

## `plant-srna-annotation/` (prototype, 2026-11-11)

```
plant-srna-annotation/
├── SKILL.md                         map of the workflow; read first
├── assets/workflow.svg              workflow figure (vector; also describes itself in text)
└── references/
    ├── glossary.md                  biology and jargon
    ├── 01_reference_genome.md
    ├── 02_library_sourcing.md       GEO/SRA selection and its pitfalls
    ├── 03_library_processing_qc.md  trimming, tag counts, QC, freezing
    ├── 04_compute_and_queue.md      environments, cores, detached queue
    ├── 05_mirna_prediction.md       ShortStack 4, miRador
    ├── 06_mirna_validation_mirscore.md
    ├── 07_phas_phasis.md
    ├── 08_structural_rna.md
    ├── 09_outputs_and_interpretation.md  what every output means
    ├── 10_release.md
    └── pitfalls.md                  symptom -> cause -> fix
```

`SKILL.md` is short and always loaded when the skill is used; each reference file is read only when its step comes up, which keeps the agent's context (and token use) small. The files contain no external links, only links to each other.

### Install for Claude Code

Copy or link the folder into your skills directory:

```bash
# personal
cp -r skills/plant-srna-annotation ~/.claude/skills/
# or for one project only
cp -r skills/plant-srna-annotation <project>/.claude/skills/
```

It then triggers on requests about plant sRNA annotation, miRNA/*PHAS* prediction, sRNA library selection, or the outputs of ShortStack, miRador, miRScore and Phasis. People can also read the files directly as a manual, starting from `SKILL.md`.

### Status

A prototype: the tomato run is still in progress and Phasis performance work is ongoing, so run times and some steps will change. Numbers quoted are preliminary rice and tomato values.
