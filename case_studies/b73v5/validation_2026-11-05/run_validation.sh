#!/usr/bin/env bash
# Record of the 2026-11-05 run on Cornell BioHPC (cbsubl02); paths are those used on that server.
# Reproduce the 2026-11-05 GFF/FASTA validation of the B73v5 curated sRNA reference v1.
# NEW validation code; not the original annotation-building code. Read-only on all inputs.
set -euo pipefail
V=/workdir/thc47/gff_validation_2026-11-05
R=/workdir/thc47/B73v5_curated_sRNA_reference_v1_2026-08-27
C=/workdir/thc47/compare_maize          # existing Bowtie1 SAMs (MaizeGDB genome), used read-only
cd "$V"

# genome is a symlink so the .fai is written here, not beside the original FASTA
ln -sfn "$C/Zm-B73-REFERENCE-NAM-5.0.fa" genome/Zm-B73-REFERENCE-NAM-5.0.fa
[ -s genome/Zm-B73-REFERENCE-NAM-5.0.fa.fai ] || samtools faidx genome/Zm-B73-REFERENCE-NAM-5.0.fa

# 1. supplied class-level helper (sequence-at-coordinate, either orientation)
for p in 21 24; do
  python3 scripts/validate_gff_product_sequences.py \
    --genome genome/Zm-B73-REFERENCE-NAM-5.0.fa \
    --gff "$R/annotation/B73v5_sRNA_reference_v1.${p}_PHAS.gff3" \
    --products "$R/fasta/${p}_PHAS.phasiRNAs.fa" \
    --feature-type exon --out-prefix results/class_level/${p}_phasiRNAs | tee logs/class_level_${p}.log
done

# 2. strict feature-ID-aware validation, all classes of the combined GFF3
python3 scripts/validate_srna_reference_strict.py \
  --release "$R" \
  --genome genome/Zm-B73-REFERENCE-NAM-5.0.fa \
  --sam21 "$C/21_phasiRNAs_Zm-B73-REFERENCE-NAM-5.0.sam" \
  --sam24 "$C/24_phasiRNAs_Zm-B73-REFERENCE-NAM-5.sam" \
  --outdir results/strict | tee logs/strict.log

{ bedtools --version; samtools --version | head -1; python3 --version;
  git -C /workdir/thc47/Phasis describe --tags --always; } > logs/software_versions.txt 2>&1
