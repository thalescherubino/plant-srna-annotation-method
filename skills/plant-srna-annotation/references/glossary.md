# Glossary: the biology and jargon you need

[Back to the skill map](../SKILL.md)

## Small RNAs in plants

- **Small RNA (sRNA):** RNA of about 18-30 nucleotides (nt). In plants most come from Dicer-like enzymes (DCLs) cutting longer double-stranded RNA. The cut products are loaded into Argonaute (AGO) proteins, which use them as guides to silence matching RNA or DNA.
- **Size classes matter.** Dicer cuts at fixed lengths, so a healthy plant library shows sharp peaks at **21 nt** (and some 22 nt) and **24 nt**. A library without these peaks is suspect (wrong trimming, degraded RNA, or not a small-RNA library at all).
- **5′ nucleotide:** AGO proteins sort sRNAs partly by their first base. 21-nt miRNAs mostly start with U (AGO1); 24-nt siRNAs mostly start with A (AGO4 family). A library's 5′U share at 21 nt and 5′A share at 24 nt is a quick biological sanity check.

## miRNAs

- **miRNA:** a 20-22-nt (rarely 24-nt) sRNA cut from a single-stranded RNA that folds back on itself into a **hairpin**. It guides cleavage or translational repression of target mRNAs.
- ***MIR* gene / precursor / hairpin:** the genomic locus and its fold-back transcript. Gene names are upper case in italics (*MIR166a*); mature products are "miR166a".
- **Mature miRNA and star (miRNA\*):** Dicer cuts a **duplex** of two strands from the hairpin stem. The strand that is kept is the mature miRNA; the partner is the star. Arms are named **5p** and **3p** by their side of the hairpin.
- **2-nt 3′ overhang:** Dicer leaves two unpaired bases at the 3′ end of each strand of the duplex. A real miRNA duplex has this structure, so miRScore checks it. It cannot be checked if the mature starts at the very first base of the supplied hairpin.
- **Precision:** of all reads mapping to a hairpin, the fraction that are the mature or star (± small shifts). Real *MIR* loci are precise; degradation produces scattered reads.
- **miRBase:** the public miRNA registry. Isoforms of one hairpin are sometimes named `.1`/`.2`. Entries can be wrong or come from another cultivar or assembly.
- **Rescue (miRScore):** when the annotated duplex fails but another duplex on the same hairpin is well supported, miRScore can propose the alternative.

## Phased siRNAs

Terms below follow the Phasis paper (Cherubino Ribeiro et al. 2026, bioRxiv) and the Phasis documents.

- **phasiRNA:** secondary siRNAs produced in phase. "Secondary" because an upstream cleavage, usually guided by a miRNA, starts their production and sets the phase register. A typical pathway: a 22-nt miRNA (or a pair of miRNAs, as for *TAS3* with miR390) cleaves a *PHAS* transcript; RNA-dependent RNA polymerase 6 (RDR6) copies the cleaved fragment into double-stranded RNA; a Dicer-like enzyme then cuts it from the cleavage end in consecutive steps: DCL4 for 21 nt; for 24 nt, DCL5 in grasses (in eudicots, which lack DCL5, another Dicer such as DCL3 is thought to act). Their 5′ ends therefore fall every 21 (or 24) nt: they are **in phase**. miRNAs are also Dicer products, but they are not phased and do not define *PHAS* loci.
- **Triggers differ by lineage.** miR482/miR2118 trigger 21-nt phasiRNAs from NB-LRR disease-resistance genes in many eudicots, including tomato; miR2118 and miR2275 trigger the large reproductive 21- and 24-nt families of grass anthers. Tomato 24-nt phasiRNAs, reported in meiotic inflorescences, have no known miRNA trigger and show unusual 12-nt phasing (Xia et al. 2019, Nature Communications). Published tomato counts are about 90-140 21-*PHAS* loci (for example 108 in Baldrich et al. 2022, Plant Physiology) and 137 24-*PHAS* loci (Xia et al. 2019).
- ***PHAS* locus:** the genomic region producing phasiRNAs; italic because it names a genomic locus, like *MIR*. It can be protein-coding (NB-LRR, PPR, MYB genes) or non-coding (*TAS* genes, long non-coding transcripts). **21-*PHAS*** loci occur in vegetative tissues and flowers, including the large reproductive 21-nt families of grass anthers; **24-*PHAS*** loci are typical of anthers and meiotic tissues. Reproductive 24-*PHAS* loci depend on the miR2275-DCL5 pathway: Arabidopsis lacks this pathway, and in maize *dcl5* mutant anthers 24-*PHAS* signal is strongly depleted. A well-supported *PHAS* locus shows the expected size class (21 or 24 nt), a coherent phased register and support at several phased positions. High sRNA abundance or apparent periodicity alone is not enough.
- **Register / phase position:** a series of expected read positions spaced by the phase length (21 or 24 nt). Reads at those positions support phasing. Reads 1 nt off are common, because Dicer processing varies a little. Phasis's main Howell score tolerates ±1 nt at the 5′ end, separately on each strand; its strict score does not. A candidate with no read at the exact register position at its anchor (HPSP) is not called *PHAS*, even if ±1-nt reads are present. Paired reads on the opposite strand sit at a fixed offset set by Dicer's 2-nt 3′ overhang. Phasis checks this geometry (an "opposite-strand partner"), but the check does not tell which strand was transcribed or targeted by the trigger.
- **Phasis score vs Howell score:** two different measures (Phasis `docs/scoring.md`). The **Phasis score** is statistical. In each window, a hypergeometric test (are reads enriched at phased positions?) and a Wilcoxon rank-sum test (are phased positions more abundant?) are combined by Stouffer's method; the windows of a locus are then combined by Fisher's method. The score is -log10 of the result, capped at 300. The **Howell score** measures how regularly phase-length reads are spaced, from their abundance and the number of occupied phase cycles in a ten-cycle window (adapted from Howell et al. 2007). Phasis reports a main (1-nt tolerant) peak, a strict (exact positions only) peak and an exact-support score. ShortStack 3 (PhaseScore) and PhaseTank (Phased_Score) use their own adaptations of the Howell score, so their values are not directly comparable with Phasis's.
- **Phased unit:** Phasis's term for a group of supported phased windows describing one local phasiRNA-producing region. The **main unit** is anchored at the **HPSP** (highest phasing score position: the register with the highest Howell score). Same-register signal connected to it is a **coherent extension**; compatible opposite-strand signal is an accepted **partner**. Secondary units are either supported windows outside the main unit (separate local peaks) or **overlapping alternative windows** (overlapping the main unit, but in another or shifted register). Phasis does not treat a secondary unit as a separate locus or precursor.
- ***PHAS*-like:** in Phasis, a candidate with detectable phased structure that lacks the register-level evidence for a high-confidence call. Typical reasons: weak exact-register support, a crowded local phasing context, competing alternative registers, or signal that is not coherent enough. Kept separately for follow-up and not added to the main calls table; analyse it apart from the *PHAS* calls. Phasis decides the class per library, so one locus can be *PHAS* in one library and *PHAS*-like in another.
- **Individual vs pooled runs:** Phasis analyses each library individually by default; pooled runs (`--pool_libraries`) combine all libraries into one virtual library first. The two modes are complementary, not ranked. Individual runs keep sample-specific signal; pooling recovers weak signal spread over many libraries but can lower precision. Pooling is not a strict superset: in the paper's maize benchmark, individual runs found 46 21-*PHAS* and 14 24-*PHAS* reference loci that the pooled run missed. In our runs, pooling found weak loci but also added dense 24-nt background.
- **hc-siRNA (heterochromatic siRNA):** abundant 24-nt siRNAs from transposons and repeats, guiding DNA methylation. They are not phased, but dense 24-nt clusters can show high abundance or apparent periodicity, the main source of false 24-*PHAS* calls.

## Structural RNAs

- **Structural ncRNAs:** tRNA, rRNA, snRNA, snoRNA. Their fragments are abundant in sRNA libraries but they are not sRNA genes; annotating them lets you separate them from real sRNA loci.
- **Covariance model (CM):** a statistical model of an RNA family's sequence plus secondary structure. **Rfam** is the library of these models; **Infernal** (`cmscan`) searches a genome with them.
- **GA (gathering) threshold:** the per-family score cut-off curated by Rfam (`--cut_ga`).
- **Clan:** a group of related Rfam families; overlapping hits from one clan compete and only the best is kept.
- **Model coverage:** how much of the model a hit spans. Full-length hits (≥ 80%) are canonical; partial hits are often repeat-derived fragments.
- **Organelle-type ncRNA / NUPT / NUMT:** bacterial-like rRNA or tRNA. Either a true organelle gene (if organelle contigs are in the assembly) or a copy of chloroplast (NUPT) or mitochondrial (NUMT) DNA inserted in the nucleus.
- **CCA tail:** tRNAs gain a CCA at the 3′ end after transcription, so mature tRNA reads may not match the genome exactly.

## Data formats and processing

- **SRA / GEO:** NCBI archives. GEO holds processed files and sample descriptions; SRA holds raw reads. An **SRX** is one sequencing experiment (usually one library), which can have several **SRR** runs.
- **Adapter:** a synthetic sequence ligated to the RNA for sequencing. sRNAs are shorter than reads, so every read contains the 3′ adapter, which must be trimmed.
- **4N / random ends, UMIs, barcodes:** some kits add random or index bases to the RNA ends. Unless removed, these make reads unmappable.
- **Untemplated 3′ additions (tailing):** cells add a few non-genomic bases (often U or A) to sRNA 3′ ends. Do not trim reads to "fix" them; that shifts the 21/24-nt peaks.
- **Tag count:** a collapsed library: each distinct sequence once, with its read count.
- **Exact mapping / multimapping:** a read matching the genome with zero mismatches; a read matching several places (common for 24-nt repeat siRNAs).
- **Bowtie 1, HISAT2:** short-read aligners used inside the tools.
- **GFF3:** the standard annotation file format (one feature per line with coordinates, strand and attributes).
