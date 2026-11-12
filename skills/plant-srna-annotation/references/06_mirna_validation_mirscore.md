# Step 6: Score miRNA candidates with miRScore

[Back to the skill map](../SKILL.md) · Previous: [05_mirna_prediction.md](05_mirna_prediction.md) · Next: [07_phas_phasis.md](07_phas_phasis.md)

## Goal

Score every candidate set **on the same libraries with the same rules**: ShortStack, miRador and the species' miRBase entries. Each hairpin gets Pass or Fail with reasons.

miRScore is not a predictor; it checks supplied candidates against plant miRNA criteria: duplex structure with a 2-nt 3′ overhang, few mismatches, mature/star read support, precision, and support in more than one library.

## Settings

- Plant mode (`-kingdom plant`), 0 mismatches, **`-rescue` always on** (it recovers real *MIR* loci whose annotated duplex is not the dominant one).
- Libraries: the frozen tag counts expanded to **uncompressed** FASTA, one record per read (Phasis `support_scripts/TagToFasta.py`). Do not include a pooled library together with the individual libraries it was made from.
- miRBase is a **tested set**, not a positive control: it contains entries that fail on real data.

## Build the inputs (`scripts/mirna/build_mirscore_inputs.py`)

```bash
build_mirscore_inputs.py --genome genome.fa --set shortstack --input shortstack_denovo/results --out shortstack/inputs
build_mirscore_inputs.py --genome genome.fa --set mirador --input mirador/output/finalAnnotatedCandidates.csv --out mirador/inputs
build_mirscore_inputs.py --genome genome.fa --set mirbase --input <miRBase dir> --species-prefix sly- --out mirbase/inputs
```

- **Predictor hairpins:** duplex span ± 20 nt, capped so no hairpin exceeds 300 nt (miRScore's plant mode fails longer precursors). Both predictors get the same rule because precision depends on hairpin span.
- **miRBase hairpins:** kept as deposited. If a mature or star sits at a hairpin end, the hairpin is located exactly in the genome and extended under the same cap.
- **Not scorable:** if such a hairpin is not found exactly in the genome (it comes from another cultivar or assembly), it cannot be extended. It is dropped with its matures and logged in `not_scorable.tsv`. Tomato: 5 of 112 (sly-MIR156b, 156c, 159, 162, 319a); BLAST showed 1-4 mismatches or a 1-nt gap against the genome.
- **miRBase isoforms:** run `hairpinHelper` on the miRBase set. It makes one hairpin entry per `.1`/`.2` mature isoform, because miRScore matches mature to hairpin by name. If no hairpin needs it (tomato), hairpinHelper writes **nothing**; then use the original hairpins.

## Run in chunks (`scripts/mirna/mirscore_chunks.py`)

miRScore runs single-threaded for much of its work. Split each set into 8 hairpin chunks, run them in parallel with 4 threads each, then merge:

```bash
mirscore_chunks.py split <set_dir> --chunks 8
# per chunk:
miRScore -fastq libs/*.fa -mature inputs/mature.fa -hairpin inputs/hairpin.fa \
  -kingdom plant -rescue -threads 4 -n <set>_<chunk> -out results
mirscore_chunks.py merge <set_dir> --name <set>
```

Chunking does not change results: each hairpin is scored independently and Bowtie reports every alignment. The merge checks that every input hairpin appears exactly once.

- Count a chunk as done **only when its results CSV exists**; miRScore can exit with code 0 and write nothing.
- If you change the input hairpin list, re-split and rerun **all** chunks. Chunk membership changes, so old chunk results no longer line up.

## Typical results (preliminary)

| Set | Rice | Tomato |
|---|---|---|
| ShortStack (discovery run) | 96/97 Pass | 106/107 Pass |
| miRador | 437/580 Pass (36 via rescue) | 343/481 Pass |
| miRBase 23.1 | 276/516 Pass (179 only via rescue) | 107 scored after 5 dropped (running) |

ShortStack candidates nearly all pass because ShortStack applies similar criteria. miRador proposes more, including many 24-nt candidates; see [09_outputs_and_interpretation.md](09_outputs_and_interpretation.md) for how to read the flags.

## When it fails

See [pitfalls.md](pitfalls.md), section miRScore: no results despite exit code 0, a chunk stopped by "miRNA starts at position 1", `.gz` libraries read as FASTQ, a slow silent rescue step.
