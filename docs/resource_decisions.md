# Adapting to incomplete resources

For every target species, search authoritative assembly and sequence sources, study metadata and the literature; record search date, versions, query, results, license and missing resources. Verify source availability at the time of use. Select the following paths according to observed evidence rather than automatically substituting maize thresholds.

| Missing or incomplete resource | Reasonable next action | Limitation to report |
| --- | --- | --- |
| Species miRNA annotation | Generate de novo candidates using miRador and ShortStack 4.1.0+ (`--dn_mirna`) with suitable small-RNA libraries and genome; evaluate compatible precursor/product models with miRScore. | Homology alone cannot establish a bona-fide MIR; without suitable libraries, leave unvalidated candidates provisional. |
| Structural RNA entries | Predict tRNAs with tRNAscan-SE; use Infernal and Rfam to search for snRNA/snoRNA/rRNA family models, with RNAcentral sequences and genome mapping where available. | Model sensitivity, genome completeness and family assignment can limit detection; report failed and ambiguous predictions. |
| Species *PHAS* catalog | Run Phasis 2.8.4 or newer for *PHAS* and *PHAS*-like discovery from quality-controlled libraries, including individual and pooled evidence where justified; review register and architecture. | A lack of suitable libraries prevents strong tissue-specific claims; pooled-only evidence is not replicate-specific evidence. |
| High-quality genome | Assess contiguity, contamination, organelle representation, ploidy and sequence completeness; consider a documented transcript-based exploratory analysis if appropriate. | Transcript-coordinate results and missing genomic loci cannot be interpreted as a complete genome annotation. |
| Independent validation dataset | Use available replication, alternate evidence and manual review; report the missing validation explicitly. | Do not label computational candidates experimentally validated. |

Tools, algorithms, model releases and databases change. An agent adapting this method should inspect current official documentation and primary sources, evaluate organism and input suitability, pin the selected versions, and record each decision. Report which branches could not be completed and why. Avoid implying that computationally predicted *PHAS* or miRNA products demonstrate mobility or biological function.
