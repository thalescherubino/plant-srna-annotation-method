# Software setup and first-use guide

Use separate environments: miRador has older dependencies, while current Phasis has its own tested Python stack. This guide gives a starting configuration, not a cross-platform lockfile; channel availability and compatibility can change. Use Linux or macOS supported by the tools, install a Conda-compatible package manager, and execute installations and analyses under your institution's computing policies. Do not install into an existing shared or production environment. A computer with sufficient storage and memory is needed for genomes, indexes and libraries.

The commands below assume you are in a clone of this methods repository. They create three named environments; check for existing names first with `conda env list` to avoid modifying someone else's environment. `conda activate` requires a configured Conda shell. On shared systems select writable environment and package-cache locations according to local policy. Record the final dependency list and platform; do not assume an environment solve that works on one system is portable to another.

## 1. Phasis: new *PHAS* discovery

Use [Phasis](https://github.com/blakemeyerslab/Phasis) **2.8.4 or newer**. The upstream 2.8.4 installation guidance uses Conda for dependencies and pip for the source checkout. The YAML follows its dependency versions; it does not install Phasis itself.

```bash
conda env create -f envs/phasis.yml
conda activate plant-srna-phasis
git clone https://github.com/blakemeyerslab/Phasis.git
cd Phasis
git rev-parse HEAD
python -m pip install -e .
phasis -version
phasis -h
python -m pip check
hisat2 --version
samtools --version
python -c "from phasis.deps_check import require_dependencies; require_dependencies()"
```

Confirm the reported version is at least 2.8.4 **before** discovery, and record the version and `git rev-parse HEAD` with the analysis. A moving `main` is not a fixed release: after verifying compatibility, preserve the commit and environment export for the study. Editable installation uses this checkout, so keep it unchanged and available. Return to the methods repository before running the remaining setup commands. Check the [upstream instructions](https://github.com/blakemeyerslab/Phasis#installation) for updated dependencies and troubleshooting.

Phasis accepts quality-controlled FASTQ, FASTA, or tag-count libraries as documented upstream; it does not substitute for adapter trimming or UMI handling. Start with individual libraries, then run pooled discovery as a separate analysis when appropriate. Select phase lengths appropriate to the organism and study; 21 and 24 are examples, not universal defaults. Use a distinct working directory per concurrent run: separate output directories alone do not isolate the shared Phasis cache. See the [Phasis examples](https://github.com/blakemeyerslab/Phasis#use-your-own-data) for exact commands and input rules.

For example, after replacing the placeholder paths and setting an appropriate CPU count, run the following **sequentially** from a new study directory with `plant-srna-phasis` active:

```bash
phasis -libs /path/to/libraries/*.tag -libformat T \
  -reference /path/to/genome.fa -phase 21 -cores 4 --outdir 21_individual
phasis -libs /path/to/libraries/*.tag -libformat T \
  -reference /path/to/genome.fa -phase 24 -cores 4 --outdir 24_individual
phasis -libs /path/to/libraries/*.tag -libformat T \
  -reference /path/to/genome.fa -phase 21 -cores 4 --pool_libraries --outdir 21_pooled
phasis -libs /path/to/libraries/*.tag -libformat T \
  -reference /path/to/genome.fa -phase 24 -cores 4 --pool_libraries --outdir 24_pooled
```

Confirm that the glob includes exactly the intended libraries; tag counts must retain the original read abundances. Inspect the classification evidence, register-supported products, per-library calls and diagnostic plots. Preserve pooled evidence separately. Review actual output and tool documentation before integrating candidates.

## 2. miRador: de novo miRNA prediction

[miRador](https://github.com/rkweku/miRador) predicts plant miRNA candidates from a genome and small-RNA data. Use the tool's own environment in a **separate shell/environment** from Phasis and miRScore:

```bash
conda deactivate
git clone --recursive https://github.com/rkweku/miRador.git
cd miRador
git rev-parse HEAD
git submodule status --recursive
conda env create -n plant-srna-mirador -f environment.yml
conda activate plant-srna-mirador
```

Follow the [miRador configuration and usage instructions](https://github.com/rkweku/miRador#configuration) to set the genome, prepared libraries, organism and output directory in a study-specific `miRador.ini`. Its documented library inputs include chopped/tag files; check format and preserve read abundance. Use a new output directory: the upstream documentation warns that existing output may be overwritten. Save the configuration, version/commit, library manifest, commands and candidate outputs. Validate the environment with the upstream test or example before a full study. If the older dependency stack cannot be solved on your platform, document the failure and use an upstream-supported installation or separately justified predictor rather than silently claiming miRador was run. Return to the methods repository for the next section.

## 3. miRScore: evaluate miRNA candidates

[miRScore](https://github.com/Aez35/miRScore) evaluates supplied miRNA candidates; it is not a de novo genome-wide predictor. The separate environment pins 0.3.7, the version documented for the B73 case:

```bash
conda deactivate
conda env create -f envs/mirscore.yml
conda activate plant-srna-mirscore
miRScore -help
bowtie --version
samtools --version
RNAfold --version
```

Prepare matched mature and precursor FASTA identifiers, with each precursor and product checked against the chosen genome. Supply separate, adapter-trimmed, quality-controlled small-RNA libraries, and select plant criteria; do not include a pooled library alongside the individual libraries that comprise it. Adapt this example to real paths and confirm flags with the installed version's help:

```bash
miRScore -mature candidates.mature.fa -hairpin candidates.precursors.fa \
  -fastq library_1.fastq.gz library_2.fastq.gz -kingdom plant
```

Inspect `miRScore_Results.csv`, per-library read evidence and flags, and any alternative models; computational pass/fail is not experimental validation. Version-specific `-rescue` alternatives require review. If no suitable libraries or matching mature/precursor inputs exist, report the candidate as unassessed or provisional. Consult the [miRScore README](https://github.com/Aez35/miRScore) for input formats and the treatment of alternative candidates.

## 4. Structural RNAs and optional comparisons

The following environment provides [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) and [Infernal](https://github.com/EddyRivasLab/infernal); family identification also requires an appropriate versioned [Rfam](https://rfam.org/) model collection and genome or candidate sequences. Record the database release, model thresholds and commands. Classify tRNAs, rRNAs, snRNAs and snoRNAs with genome evidence and explicitly documented inclusion rules.

```bash
conda deactivate
conda env create -f envs/structural.yml
conda activate plant-srna-structural
tRNAscan-SE -h
cmscan -h
```

[ShortStack](https://github.com/MikeAxtell/ShortStack) is optional for complementary *PHAS* evidence. Install it separately following the documentation for the chosen version and record its mapping and calling parameters. Its results are not interchangeable with Phasis results without an input, assembly and method audit.

## 5. Freeze provenance

After verifying each tool, save exact versions, source commit IDs where applicable, `conda list --explicit` for each environment, database releases, input checksums, configuration files, commands and logs. `conda list --explicit` is platform-specific and does not capture external source checkouts. Run representative upstream examples on the target platform before interpreting real datasets. The environments above are starting specifications and have not been tested by this repository on every platform.
