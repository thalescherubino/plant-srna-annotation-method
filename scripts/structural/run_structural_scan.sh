#!/usr/bin/env bash
# Genome-wide structural RNA prediction, applied identically to every species (policy decided 2026-11-06):
#   1. tRNAscan-SE (eukaryotic mode) + EukHighConfidenceFilter (flags, does not remove);
#   2. Infernal cmscan with the Rfam structural subset (rRNA, snRNA, snoRNA, scaRNA), following the Rfam
#      genome-annotation recipe: --cut_ga --rfam --nohmmonly --fmt 2 --clanin, Z = 2 x genome size in Mb.
# RNAcentral records are compared AFTER this prediction, not used as its starting point.
#
# Usage: run_structural_scan.sh <genome.fa (faidx-indexed)> <rfam_dir> <outdir> <prefix> [threads]
# Requires the plant-srna-structural environment (envs/structural.yml) to be active, and
# scripts/structural/make_rfam_structural_subset.py to have been run on <rfam_dir>.
set -euo pipefail
G=$1; RF=$2; O=$3; P=$4; T=${5:-8}
[ -s "$G.fai" ] || { echo "missing $G.fai (run samtools faidx)" >&2; exit 1; }
[ -s "$RF/Rfam.structural_subset.cm.i1m" ] || { echo "missing pressed structural subset in $RF" >&2; exit 1; }
mkdir -p "$O/trnascan/hc" "$O/rfam_cmscan"
Z=$(awk '{t+=$2} END{printf "%.6f", 2*t/1e6}' "$G.fai")
echo "genome=$G Z=$Z threads=$T" | tee "$O/scan_parameters.txt"
tRNAscan-SE -h 2>&1 | sed -n 2p >> "$O/scan_parameters.txt"
cmscan -h | sed -n 2p >> "$O/scan_parameters.txt"

( cd "$O/trnascan" && tRNAscan-SE -E -H --detail --thread "$T" -o "$P.trnascan.out" -f "$P.trnascan.ss" \
    -s "$P.trnascan.iso" -m "$P.trnascan.stats" -b "$P.trnascan.bed" -j "$P.trnascan.gff3" -a "$P.trnascan.fa" \
    -l "$P.trnascan.log" "$G" && \
  EukHighConfidenceFilter -i "$P.trnascan.out" -s "$P.trnascan.ss" -o hc -p "$P.trnascan.hc" )

( cd "$O/rfam_cmscan" && cmscan -Z "$Z" --cut_ga --rfam --nohmmonly --cpu "$T" --fmt 2 \
    --clanin "$RF/Rfam.structural_subset.clanin" --tblout "$P.rfam_structural.tblout" \
    -o "$P.rfam_structural.cmscan" "$RF/Rfam.structural_subset.cm" "$G" )
echo "done $(date -Is)" >> "$O/scan_parameters.txt"
