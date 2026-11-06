# Case study: maize B73 RefGen_v5

The documented reference used `Zm-B73-REFERENCE-NAM-5.0`. The assembly FASTA SHA-256 recorded in the source report was `a40ee35888a9c0b956c448752f59637639a7906d7fa14971c9047eb088ded2c2`. The released genome annotation and sequence references were constructed from distinct evidence streams and reconciled before packaging.

| Branch | Recorded approach |
| --- | --- |
| miRNA | miRBase v23, the miRScore study's supplementary S4 de novo candidates from miRador, ShortStack and miRDeep-P2, and expanded MIR2275 candidates; plant-mode miRScore 0.3.7, commit `64c3abfe93157a849506fe49c78cddb401e7ac05`; exact and repaired precursor models reviewed and passing models merged by same-strand precursor overlap. S4 was a published source; this does not establish that all three predictors were rerun in this curation. |
| Structural RNAs | tRNAscan-SE 2.0.13 in eukaryotic mode predicted genome-encoded tRNAs; RNAcentral v27 provided snRNA/snoRNA/rRNA candidates, Bowtie 1.3.1 supplied exact full-length genomic mappings, and Infernal 1.1.4 `cmscan --cut_ga --rfam` classified candidates with Rfam covariance models. tRNA pseudogenes were excluded; separate genomic and CCA-aware tRNA references were retained. |
| *PHAS* | 49 quality-controlled 18–30-nt tag libraries; Phasis 2.8.3 commit `00fb8c0e6cf9c332a2720cf0d31f540cccc35331`, individual and pooled 21/24-nt discovery, with ShortStack 3.8.5 comparison and manual review. Exact mapping, maxhits 25 and minimum Howell score 12.5 were case settings, not general defaults. |
| Reconciliation | MIR precedence over overlapping *PHAS*; specific exact structural-RNA containment exclusions; maximum Peak Howell score for overlapping cross-phase calls. All are release-specific decisions preserved in an exclusion ledger. |

The released resource reported 220 bona-fide MIR loci, 439 21-*PHAS* and 452 24-*PHAS* loci, plus structural RNA classes and explicitly exploratory candidate sets. Some accepted *PHAS* additions were locus-only; no phasiRNA products were invented for them. Counts reflect this particular release and are not expectations for another genome. The historical 2.8.3 version is retained here for provenance; new discovery in this method requires Phasis 2.8.4 or newer.

## Historical implementation

The [scripts](scripts/) are copies of five local case-specific builders used or associated with curation and release. The unified GFF3 builder is directly named in the release manifest. For the other scripts, association with this release needs individual command/log verification; their presence alone is not proof of execution. The source release manifest records input checksums but does not establish a build-time checksum for every script. These copies reflect the recovered local versions, not verified immutable executed binaries. They assume B73-specific inputs, coordinate conventions and decisions; they are not a ready-to-run pipeline for other species.

## Sources

- B73 assembly: [Hufford et al. 2021](https://doi.org/10.1126/science.abg5289).
- miRScore: [Vanek et al. 2025](https://doi.org/10.1371/journal.pcbi.1013663).
- [RNAcentral](https://doi.org/10.1093/nar/gku991), [Rfam](https://doi.org/10.1093/nar/gkaa1047), [tRNAscan-SE](https://doi.org/10.1093/nar/gkab688).
- [Phasis](https://doi.org/10.1101/158832); [ShortStack](https://doi.org/10.1534/g3.116.030452).
