#!/usr/bin/env python3
"""Build exploratory sRNA classes and audit cross-class conflicts in B73 v5.

This script does not change the strict bona-fide MIR or PHAS releases.  It:

* creates a sequence-level ``miRNA_unsupported_candidate`` reference from
  source-authorized 19--24-nt products of candidates excluded from the strict
  miRScore-validated MIR release;
* creates a manually delimited ``24_PHAS_like`` class from the 113 pooled-only
  loci rejected during plot review, without treating their products as
  validated phasiRNAs;
* audits genomic interval overlaps among MIR, PHAS, PHAS-like, tRNA, rRNA,
  snRNA, and snoRNA features; and
* audits exact short-product conflicts and exact containment of short products
  in structural-ncRNA reference sequences on either strand.

The sequence-equivalence FASTAs are mapping-safe class references.  Exact
sequences already represented by a stricter class are retained in ledgers but
are not duplicated in a lower-confidence sequence-equivalence FASTA.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shlex
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


VERSION = "1.0.0"
DNA_RE = re.compile(r"^[ACGT]+$")
SHORT_MIN = 19
SHORT_MAX = 24


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mirna-candidate-ledger", type=Path, required=True)
    parser.add_argument("--mirna-candidate-arms", type=Path, required=True)
    parser.add_argument("--mirna-coordinate-proposals", type=Path, required=True)
    parser.add_argument("--mirna-local-alignments", type=Path, required=True)
    parser.add_argument("--mirna-bona-fide-loci", type=Path, required=True)
    parser.add_argument("--mirna-bona-fide-products", type=Path, required=True)
    parser.add_argument("--phas-manual-decisions", type=Path, required=True)
    parser.add_argument("--phas21-loci", type=Path, required=True)
    parser.add_argument("--phas24-loci", type=Path, required=True)
    parser.add_argument("--phas21-products", type=Path, required=True)
    parser.add_argument("--phas24-products", type=Path, required=True)
    parser.add_argument("--structural-membership", type=Path, required=True)
    parser.add_argument("--trna-cca-membership", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    parser.add_argument(
        "--release-label", default="B73v5_sRNA_candidate_classes_v1_2026-08-26"
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dna(sequence: str) -> str:
    return sequence.upper().replace("U", "T")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def natural_contig_key(contig: str) -> tuple[int, int | str]:
    match = re.fullmatch(r"chr(\d+)", contig, flags=re.IGNORECASE)
    return (0, int(match.group(1))) if match else (1, contig)


def normalize_contig(contig: str) -> str:
    return f"chr{contig}" if contig.isdigit() else contig


def split_values(value: str, separator: str = ",") -> list[str]:
    return [item for item in value.split(separator) if item]


def write_fasta(path: Path, records: Iterable[tuple[str, str, str]]) -> int:
    count = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8") as handle:
        for identifier, sequence, description in records:
            handle.write(f">{identifier}")
            if description:
                handle.write(f" {description}")
            handle.write("\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")
            count += 1
    return count


@dataclass(frozen=True)
class Interval:
    feature_class: str
    feature_id: str
    contig: str
    start: int
    end: int
    strand: str
    evidence: str


def parse_structural_loci(raw: str, sequence_length: int):
    """Yield B73-v5 placements encoded as contig:start[:strand] or ranges."""
    for item in raw.split(";"):
        if not item:
            continue
        match = re.fullmatch(r"([^:]+):(\d+)(?:-(\d+))?:([+\-])", item)
        if not match:
            raise ValueError(f"Unrecognized exact_B73_loci item: {item!r}")
        contig, start, explicit_end, strand = match.groups()
        start_i = int(start)
        end_i = int(explicit_end) if explicit_end else start_i + sequence_length - 1
        yield normalize_contig(contig), start_i, end_i, strand


def build_unsupported_mirna(
    args: argparse.Namespace,
    supported_sequences: set[str],
) -> tuple[list[dict[str, object]], list[Interval], dict[str, set[str]], dict[str, int]]:
    ledger = {row["internal_id"]: row for row in read_tsv(args.mirna_candidate_ledger)}
    excluded = {
        candidate_id
        for candidate_id, row in ledger.items()
        if row["final_release_status"] == "excluded"
    }
    arm_rows = read_tsv(args.mirna_candidate_arms)
    products: list[dict[str, object]] = []
    invalid_counts: Counter[str] = Counter()
    seen_products: set[tuple[str, str, str]] = set()
    for row in arm_rows:
        candidate_id = row["internal_id"]
        if candidate_id not in excluded:
            continue
        sequence = dna(row["mature_sequence"])
        if not DNA_RE.fullmatch(sequence):
            invalid_counts["non_ACGT_product"] += 1
            continue
        if not SHORT_MIN <= len(sequence) <= SHORT_MAX:
            invalid_counts["outside_19_24_nt"] += 1
            continue
        key = (candidate_id, row["arm"], sequence)
        if key in seen_products:
            continue
        seen_products.add(key)
        source = ledger[candidate_id]
        canonical_conflict = sequence in supported_sequences
        products.append(
            {
                "unsupported_product_id": "",
                "internal_id": candidate_id,
                "display_name": source["display_name"],
                "family": source["family"],
                "arm": row["arm"],
                "mature_sequence_DNA": sequence,
                "length_nt": len(sequence),
                "source": source["source"],
                "original_id": source["original_id"],
                "predictor": source["predictor"],
                "exclusion_reason": source["exclusion_reason"],
                "source_role": row.get("source_role", ""),
                "identical_to_bona_fide_miRNA_product": "yes" if canonical_conflict else "no",
                "sequence_reference_status": (
                    "ledger_only_identical_to_bona_fide_miRNA"
                    if canonical_conflict
                    else "include_in_unsupported_candidate_sequence_reference"
                ),
            }
        )
    products.sort(key=lambda row: (str(row["internal_id"]), str(row["arm"]), str(row["mature_sequence_DNA"])))
    for index, row in enumerate(products, start=1):
        row["unsupported_product_id"] = f"ZmaB73v5-miRUC-P{index:06d}"

    all_product_path = args.output_dir / f"{args.release_label}.miRNA_unsupported_candidate.products.tsv"
    product_fields = [
        "unsupported_product_id",
        "internal_id",
        "display_name",
        "family",
        "arm",
        "mature_sequence_DNA",
        "length_nt",
        "source",
        "original_id",
        "predictor",
        "exclusion_reason",
        "source_role",
        "identical_to_bona_fide_miRNA_product",
        "sequence_reference_status",
    ]
    write_tsv(all_product_path, products, product_fields)

    by_sequence: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in products:
        if row["sequence_reference_status"] == "include_in_unsupported_candidate_sequence_reference":
            by_sequence[str(row["mature_sequence_DNA"])].append(row)
    equivalence_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str, str]] = []
    for index, sequence in enumerate(sorted(by_sequence), start=1):
        members = by_sequence[sequence]
        sequence_id = f"ZmaB73v5-miRUC-SEQ{index:06d}"
        candidate_ids = sorted({str(row["internal_id"]) for row in members})
        equivalence_rows.append(
            {
                "sequence_reference_id": sequence_id,
                "mature_sequence_DNA": sequence,
                "length_nt": len(sequence),
                "candidate_count": len(candidate_ids),
                "candidate_ids": ",".join(candidate_ids),
                "display_names": ",".join(sorted({str(row["display_name"]) for row in members})),
                "families": ",".join(sorted({str(row["family"]) for row in members if row["family"]})),
                "sources": ",".join(sorted({str(row["source"]) for row in members})),
                "status": "unsupported_MIR_candidate_sequence_not_miRScore_validated",
            }
        )
        fasta_records.append(
            (
                sequence_id,
                sequence,
                f"class=miRNA_unsupported_candidate candidate_count={len(candidate_ids)}",
            )
        )
    equivalence_path = args.output_dir / f"{args.release_label}.miRNA_unsupported_candidate.sequence_equivalence.tsv"
    equivalence_fields = [
        "sequence_reference_id",
        "mature_sequence_DNA",
        "length_nt",
        "candidate_count",
        "candidate_ids",
        "display_names",
        "families",
        "sources",
        "status",
    ]
    write_tsv(equivalence_path, equivalence_rows, equivalence_fields)
    write_fasta(
        args.output_dir / f"{args.release_label}.miRNA_unsupported_candidate.sequence_equivalence.fa",
        fasta_records,
    )

    coordinate_rows: list[dict[str, object]] = []
    intervals: list[Interval] = []
    coordinate_candidates: set[str] = set()
    for row in read_tsv(args.mirna_coordinate_proposals):
        candidate_id = row["internal_id"]
        if candidate_id not in excluded or row["coordinate_acceptance"] != "accepted":
            continue
        coordinate_id = f"miRUC_{row['projection_id']}"
        contig = normalize_contig(row["chromosome"])
        start = int(row["start_1based"])
        end = int(row["end_1based"])
        coordinate_candidates.add(candidate_id)
        coordinate_rows.append(
            {
                "coordinate_model_id": coordinate_id,
                "internal_id": candidate_id,
                "display_name": row["display_name"],
                "family": row["family"],
                "chromosome": contig,
                "start_1based": start,
                "end_1based": end,
                "strand": row["strand"],
                "coordinate_method": row["coordinate_method"],
                "coordinate_status": "B73v5_compatible_but_not_miRScore_release_qualified",
                "exclusion_reason": ledger[candidate_id]["exclusion_reason"],
            }
        )
        intervals.append(
            Interval(
                "miRNA_unsupported_candidate",
                coordinate_id,
                contig,
                start,
                end,
                row["strand"],
                row["coordinate_method"],
            )
        )

    # Add confident local placements only for excluded input-incompatible records
    # that did not already receive an accepted coordinate proposal.
    for row in read_tsv(args.mirna_local_alignments):
        candidate_id = row["internal_id"]
        if (
            candidate_id not in excluded
            or candidate_id in coordinate_candidates
            or not row["local_alignment_status"].startswith("confident_full")
        ):
            continue
        for placement_number, location in enumerate(split_values(row["best_locations"], ";"), start=1):
            match = re.fullmatch(r"([^:]+):(\d+)-(\d+):([+\-])", location)
            if not match:
                raise ValueError(f"Unrecognized miRNA local placement: {location!r}")
            contig, start, end, strand = match.groups()
            coordinate_id = f"miRUC_LOCAL_{candidate_id}_{placement_number}"
            coordinate_rows.append(
                {
                    "coordinate_model_id": coordinate_id,
                    "internal_id": candidate_id,
                    "display_name": row["display_name"],
                    "family": row["family"],
                    "chromosome": normalize_contig(contig),
                    "start_1based": int(start),
                    "end_1based": int(end),
                    "strand": strand,
                    "coordinate_method": "confident_full_query_local_alignment",
                    "coordinate_status": "review_coordinate_not_miRScore_release_qualified",
                    "exclusion_reason": ledger[candidate_id]["exclusion_reason"],
                }
            )
            intervals.append(
                Interval(
                    "miRNA_unsupported_candidate",
                    coordinate_id,
                    normalize_contig(contig),
                    int(start),
                    int(end),
                    strand,
                    "confident_full_query_local_alignment",
                )
            )
    coordinate_rows.sort(
        key=lambda row: (
            natural_contig_key(str(row["chromosome"])),
            int(row["start_1based"]),
            int(row["end_1based"]),
            str(row["coordinate_model_id"]),
        )
    )
    coordinate_path = args.output_dir / f"{args.release_label}.miRNA_unsupported_candidate.coordinate_models.tsv"
    coordinate_fields = [
        "coordinate_model_id",
        "internal_id",
        "display_name",
        "family",
        "chromosome",
        "start_1based",
        "end_1based",
        "strand",
        "coordinate_method",
        "coordinate_status",
        "exclusion_reason",
    ]
    write_tsv(coordinate_path, coordinate_rows, coordinate_fields)

    sequence_members: dict[str, set[str]] = defaultdict(set)
    for row in products:
        sequence_members[str(row["mature_sequence_DNA"])].add(str(row["unsupported_product_id"]))
    summary = {
        "excluded_source_candidates": len(excluded),
        "valid_unsupported_product_records": len(products),
        "candidates_with_valid_products": len({str(row["internal_id"]) for row in products}),
        "products_identical_to_bona_fide_miRNA": sum(
            row["identical_to_bona_fide_miRNA_product"] == "yes" for row in products
        ),
        "mapping_safe_unique_sequences": len(equivalence_rows),
        "coordinate_models": len(coordinate_rows),
        "coordinate_resolved_candidates": len({str(row["internal_id"]) for row in coordinate_rows}),
        **{f"invalid_{key}": value for key, value in invalid_counts.items()},
    }
    return products, intervals, sequence_members, summary


def build_phas_like(
    args: argparse.Namespace,
    bona_fide_phas_sequences: set[str],
) -> tuple[list[dict[str, object]], list[Interval], dict[str, set[str]], dict[str, int]]:
    decisions = [
        row
        for row in read_tsv(args.phas_manual_decisions)
        if row["manual_curation_class"] == "PHAS_like"
    ]
    if any(row["phase_nt"] != "24" for row in decisions):
        raise ValueError("The manually delimited PHAS-like cohort unexpectedly includes non-24-nt loci")
    locus_rows: list[dict[str, object]] = []
    intervals: list[Interval] = []
    pooled_to_locus: dict[str, str] = {}
    for index, row in enumerate(
        sorted(
            decisions,
            key=lambda item: (
                natural_contig_key(normalize_contig(item["chromosome"])),
                int(item["start_1based"]),
                int(item["end_1based"]),
            ),
        ),
        start=1,
    ):
        locus_id = f"ZmaB73v5-24PHASLIKE-{index:06d}"
        pooled_ids = split_values(row["pooled_identifiers"])
        if not pooled_ids:
            raise ValueError(f"PHAS-like decision lacks pooled identifier: {row['combined_identifier']}")
        for pooled_id in pooled_ids:
            if pooled_id in pooled_to_locus:
                raise ValueError(f"Duplicate PHAS-like pooled identifier: {pooled_id}")
            pooled_to_locus[pooled_id] = locus_id
        contig = normalize_contig(row["chromosome"])
        start = int(row["start_1based"])
        end = int(row["end_1based"])
        locus_rows.append(
            {
                "phas_like_locus_id": locus_id,
                "phase_nt": 24,
                "chromosome": contig,
                "start_1based": start,
                "end_1based": end,
                "source_identifier": row["combined_identifier"],
                "pooled_identifiers": row["pooled_identifiers"],
                "manual_curation_class": "PHAS_like",
                "include_in_bona_fide_PHAS_catalog": "no",
                "exploratory_class_status": "include_as_separate_24_PHAS_like_class",
                "manual_curation_basis": row["manual_curation_basis"],
                "product_semantics": "Phasis_reported_product_candidate_not_validated_phasiRNA",
            }
        )
        intervals.append(
            Interval(
                "24-PHAS_like",
                locus_id,
                contig,
                start,
                end,
                ".",
                "manual_pooled_only_plot_exclusion_from_bona_fide_PHAS",
            )
        )

    write_tsv(
        args.output_dir / f"{args.release_label}.24_PHAS_like.loci.tsv",
        locus_rows,
        [
            "phas_like_locus_id",
            "phase_nt",
            "chromosome",
            "start_1based",
            "end_1based",
            "source_identifier",
            "pooled_identifiers",
            "manual_curation_class",
            "include_in_bona_fide_PHAS_catalog",
            "exploratory_class_status",
            "manual_curation_basis",
            "product_semantics",
        ],
    )

    aggregated: dict[tuple[str, str], dict[str, object]] = {}
    raw_match_count = 0
    for row in read_tsv(args.phas24_products):
        source_identifier = row["identifier"]
        if source_identifier not in pooled_to_locus:
            continue
        raw_match_count += 1
        locus_id = pooled_to_locus[source_identifier]
        sequence = dna(row["tag_seq"])
        if not DNA_RE.fullmatch(sequence) or len(sequence) != 24:
            raise ValueError(f"Invalid 24-PHAS-like product sequence: {sequence!r}")
        key = (locus_id, sequence)
        record = aggregated.setdefault(
            key,
            {
                "phas_like_product_id": "",
                "phas_like_locus_id": locus_id,
                "mature_sequence_DNA": sequence,
                "length_nt": len(sequence),
                "source_row_count": 0,
                "reported_abundance_sum": 0.0,
                "reported_abundance_max": 0.0,
                "source_libraries": set(),
                "register_classes": set(),
                "source_identifiers": set(),
                "identical_to_bona_fide_PHAS_product": (
                    "yes" if sequence in bona_fide_phas_sequences else "no"
                ),
                "product_status": "PHAS_like_product_candidate_not_validated_phasiRNA",
            },
        )
        abundance = float(row["abun"])
        record["source_row_count"] = int(record["source_row_count"]) + 1
        record["reported_abundance_sum"] = float(record["reported_abundance_sum"]) + abundance
        record["reported_abundance_max"] = max(float(record["reported_abundance_max"]), abundance)
        record["source_libraries"].add(row["alib"])
        record["register_classes"].add(row["register_class"])
        record["source_identifiers"].add(source_identifier)
    if set(pooled_to_locus) != {
        source_id
        for record in aggregated.values()
        for source_id in record["source_identifiers"]
    }:
        missing = sorted(
            set(pooled_to_locus)
            - {
                source_id
                for record in aggregated.values()
                for source_id in record["source_identifiers"]
            }
        )
        raise ValueError(f"PHAS-like loci lack product records: {missing[:5]}")

    product_rows: list[dict[str, object]] = []
    for index, key in enumerate(sorted(aggregated), start=1):
        record = aggregated[key]
        record["phas_like_product_id"] = f"ZmaB73v5-24PHASLIKE-P{index:07d}"
        record["source_libraries"] = ",".join(sorted(record["source_libraries"]))
        record["register_classes"] = ",".join(sorted(record["register_classes"]))
        record["source_identifiers"] = ",".join(sorted(record["source_identifiers"]))
        record["reported_abundance_sum"] = f"{float(record['reported_abundance_sum']):.6f}"
        record["reported_abundance_max"] = f"{float(record['reported_abundance_max']):.6f}"
        product_rows.append(record)
    product_fields = [
        "phas_like_product_id",
        "phas_like_locus_id",
        "mature_sequence_DNA",
        "length_nt",
        "source_row_count",
        "reported_abundance_sum",
        "reported_abundance_max",
        "source_libraries",
        "register_classes",
        "source_identifiers",
        "identical_to_bona_fide_PHAS_product",
        "product_status",
    ]
    write_tsv(
        args.output_dir / f"{args.release_label}.24_PHAS_like.product_candidates.tsv",
        product_rows,
        product_fields,
    )

    by_sequence: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in product_rows:
        if row["identical_to_bona_fide_PHAS_product"] == "no":
            by_sequence[str(row["mature_sequence_DNA"])].append(row)
    equivalence_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str, str]] = []
    for index, sequence in enumerate(sorted(by_sequence), start=1):
        members = by_sequence[sequence]
        sequence_id = f"ZmaB73v5-24PHASLIKE-SEQ{index:07d}"
        loci = sorted({str(row["phas_like_locus_id"]) for row in members})
        equivalence_rows.append(
            {
                "sequence_reference_id": sequence_id,
                "mature_sequence_DNA": sequence,
                "length_nt": len(sequence),
                "locus_count": len(loci),
                "phas_like_locus_ids": ",".join(loci),
                "status": "exploratory_24_PHAS_like_sequence_not_bona_fide_phasiRNA",
            }
        )
        fasta_records.append(
            (
                sequence_id,
                sequence,
                f"class=24_PHAS_like locus_count={len(loci)}",
            )
        )
    write_tsv(
        args.output_dir / f"{args.release_label}.24_PHAS_like.sequence_equivalence.tsv",
        equivalence_rows,
        [
            "sequence_reference_id",
            "mature_sequence_DNA",
            "length_nt",
            "locus_count",
            "phas_like_locus_ids",
            "status",
        ],
    )
    write_fasta(
        args.output_dir / f"{args.release_label}.24_PHAS_like.sequence_equivalence.fa",
        fasta_records,
    )
    sequence_members: dict[str, set[str]] = defaultdict(set)
    for row in product_rows:
        sequence_members[str(row["mature_sequence_DNA"])].add(str(row["phas_like_product_id"]))
    summary = {
        "manual_24_PHAS_like_loci": len(locus_rows),
        "source_product_rows": raw_match_count,
        "locus_sequence_product_candidates": len(product_rows),
        "unique_product_candidate_sequences": len(sequence_members),
        "products_identical_to_bona_fide_PHAS": sum(
            row["identical_to_bona_fide_PHAS_product"] == "yes" for row in product_rows
        ),
        "mapping_safe_unique_sequences": len(equivalence_rows),
    }
    return product_rows, intervals, sequence_members, summary


def load_bona_fide_mirna(
    args: argparse.Namespace,
) -> tuple[list[Interval], dict[str, set[str]], set[str]]:
    intervals = [
        Interval(
            "miRNA_bona_fide",
            row["release_locus_id"],
            normalize_contig(row["chromosome"]),
            int(row["start_1based"]),
            int(row["end_1based"]),
            row["strand"],
            "strict_miRScore_validated_B73v5_MIR_v1",
        )
        for row in read_tsv(args.mirna_bona_fide_loci)
    ]
    sequence_members: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(args.mirna_bona_fide_products):
        sequence_members[dna(row["mature_sequence_DNA"])].add(row["product_id"])
    return intervals, sequence_members, set(sequence_members)


def load_bona_fide_phas(
    locus_path: Path,
    product_path: Path,
    phase: int,
) -> tuple[list[Interval], dict[str, set[str]], set[str]]:
    locus_rows = read_tsv(locus_path)
    source_to_release = {
        row["source_identifier"]: row["reference_locus_id"]
        for row in locus_rows
        if row["source_method"] == "Phasis"
    }
    intervals = [
        Interval(
            f"{phase}-PHAS_bona_fide",
            row["reference_locus_id"],
            normalize_contig(row["chromosome"]),
            int(row["start_1based"]),
            int(row["end_1based"]),
            ".",
            f"manual_bona_fide_{phase}PHAS_locus_reference",
        )
        for row in locus_rows
    ]
    sequence_members: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(product_path):
        if row["identifier"] not in source_to_release:
            continue
        sequence = dna(row["tag_seq"])
        if not DNA_RE.fullmatch(sequence) or len(sequence) != phase:
            raise ValueError(f"Invalid {phase}-PHAS product: {sequence!r}")
        sequence_members[sequence].add(source_to_release[row["identifier"]])
    return intervals, sequence_members, set(sequence_members)


def load_structural(
    args: argparse.Namespace,
) -> tuple[list[Interval], dict[str, dict[str, str]], dict[str, int]]:
    membership = read_tsv(args.structural_membership)
    intervals: list[Interval] = []
    structural_sequences: dict[str, dict[str, str]] = {}
    seen_intervals: set[tuple[str, str, str, int, int, str]] = set()
    for row in membership:
        rna_class = row["class"]
        if rna_class not in {"tRNA", "rRNA", "snRNA", "snoRNA"}:
            continue
        sequence = dna(row["canonical_dna_sequence"])
        if not DNA_RE.fullmatch(sequence):
            raise ValueError(f"Ambiguous structural sequence: {row['reference_id']}")
        structural_sequences[row["reference_id"]] = {
            "class": rna_class,
            "sequence": sequence,
        }
        for contig, start, end, strand in parse_structural_loci(
            row["exact_B73_loci"], int(row["sequence_length"])
        ):
            key = (rna_class, row["reference_id"], contig, start, end, strand)
            if key in seen_intervals:
                continue
            seen_intervals.add(key)
            intervals.append(
                Interval(
                    rna_class,
                    row["reference_id"],
                    contig,
                    start,
                    end,
                    strand,
                    "exact_full_length_B73v5_placement",
                )
            )

    # Replace baseline tRNA sequences by the selected CCA-aware sensitivity
    # reference while retaining genomic placements from the baseline IDs.
    for row in read_tsv(args.trna_cca_membership):
        sequence = dna(row["alternative_sequence"])
        if not DNA_RE.fullmatch(sequence):
            raise ValueError(f"Ambiguous CCA-aware tRNA sequence: {row['alternative_reference_id']}")
        structural_sequences[row["alternative_reference_id"]] = {
            "class": "tRNA",
            "sequence": sequence,
        }
    for reference_id in [key for key, value in structural_sequences.items() if value["class"] == "tRNA"]:
        if reference_id.startswith("zma_B73v5_tRNA_") and not reference_id.startswith("zma_B73v5_tRNA_CCA_"):
            del structural_sequences[reference_id]

    summary = Counter(interval.feature_class for interval in intervals)
    summary.update(
        {f"{rna_class}_unique_reference_sequences": sum(
            value["class"] == rna_class for value in structural_sequences.values()
        ) for rna_class in ("tRNA", "rRNA", "snRNA", "snoRNA")}
    )
    return intervals, structural_sequences, dict(summary)


def audit_genomic_overlaps(intervals: list[Interval]) -> list[dict[str, object]]:
    by_contig: dict[str, list[Interval]] = defaultdict(list)
    for interval in set(intervals):
        by_contig[interval.contig].append(interval)
    rows: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for contig in sorted(by_contig, key=natural_contig_key):
        ordered = sorted(by_contig[contig], key=lambda item: (item.start, item.end, item.feature_class, item.feature_id))
        active: list[Interval] = []
        for current in ordered:
            active = [item for item in active if item.end >= current.start]
            for other in active:
                if other.feature_class == current.feature_class:
                    continue
                first, second = sorted(
                    (other, current), key=lambda item: (item.feature_class, item.feature_id, item.start, item.end)
                )
                overlap_start = max(first.start, second.start)
                overlap_end = min(first.end, second.end)
                key = (
                    first.feature_class,
                    first.feature_id,
                    second.feature_class,
                    second.feature_id,
                    contig,
                    overlap_start,
                    overlap_end,
                )
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "class_a": first.feature_class,
                        "feature_a": first.feature_id,
                        "start_a": first.start,
                        "end_a": first.end,
                        "strand_a": first.strand,
                        "class_b": second.feature_class,
                        "feature_b": second.feature_id,
                        "start_b": second.start,
                        "end_b": second.end,
                        "strand_b": second.strand,
                        "chromosome": contig,
                        "overlap_start_1based": overlap_start,
                        "overlap_end_1based": overlap_end,
                        "overlap_nt": overlap_end - overlap_start + 1,
                        "strand_relationship": (
                            "unknown"
                            if "." in {first.strand, second.strand}
                            else "same"
                            if first.strand == second.strand
                            else "opposite"
                        ),
                        "interpretation": "review_cross_class_annotation;not_automatic_exclusion",
                    }
                )
            active.append(current)
    rows.sort(
        key=lambda row: (
            natural_contig_key(str(row["chromosome"])),
            int(row["overlap_start_1based"]),
            str(row["class_a"]),
            str(row["class_b"]),
            str(row["feature_a"]),
            str(row["feature_b"]),
        )
    )
    return rows


def audit_short_exact_conflicts(
    class_sequences: dict[str, dict[str, set[str]]],
) -> list[dict[str, object]]:
    by_sequence: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for feature_class, sequences in class_sequences.items():
        for sequence, member_ids in sequences.items():
            by_sequence[sequence][feature_class].update(member_ids)
    rows: list[dict[str, object]] = []
    for sequence in sorted(by_sequence):
        classes = by_sequence[sequence]
        if len(classes) < 2:
            continue
        rows.append(
            {
                "mature_sequence_DNA": sequence,
                "length_nt": len(sequence),
                "class_count": len(classes),
                "classes": ",".join(sorted(classes)),
                "class_memberships": ";".join(
                    f"{feature_class}={','.join(sorted(classes[feature_class]))}"
                    for feature_class in sorted(classes)
                ),
                "mapping_interpretation": "same_query_sequence_is_not_class_identifiable",
            }
        )
    return rows


def audit_structural_containment(
    class_sequences: dict[str, dict[str, set[str]]],
    structural_sequences: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    short_lookup: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for feature_class, sequences in class_sequences.items():
        for sequence, member_ids in sequences.items():
            short_lookup[sequence][feature_class].update(member_ids)
    by_length: dict[int, set[str]] = defaultdict(set)
    for sequence in short_lookup:
        by_length[len(sequence)].add(sequence)

    matches: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for reference_id, record in structural_sequences.items():
        structural_class = record["class"]
        sequence = record["sequence"]
        for length, queries in by_length.items():
            if length > len(sequence):
                continue
            for offset in range(0, len(sequence) - length + 1):
                kmer = sequence[offset : offset + length]
                if kmer in queries:
                    matches[(kmer, structural_class, reference_id)].add("forward")
                reverse = reverse_complement(kmer)
                if reverse in queries:
                    matches[(reverse, structural_class, reference_id)].add("reverse_complement")
    rows: list[dict[str, object]] = []
    for (query, structural_class, reference_id), orientations in sorted(matches.items()):
        for short_class in sorted(short_lookup[query]):
            rows.append(
                {
                    "mature_sequence_DNA": query,
                    "length_nt": len(query),
                    "short_RNA_class": short_class,
                    "short_RNA_member_ids": ",".join(sorted(short_lookup[query][short_class])),
                    "structural_class": structural_class,
                    "structural_reference_id": reference_id,
                    "orientation_in_structural_reference": ",".join(sorted(orientations)),
                    "mapping_interpretation": "exact_short_read_maps_to_structural_reference",
                }
            )
    return rows


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    mir_intervals, mir_sequences, supported_mir_sequences = load_bona_fide_mirna(args)
    phas21_intervals, phas21_sequences, phas21_unique = load_bona_fide_phas(
        args.phas21_loci, args.phas21_products, 21
    )
    phas24_intervals, phas24_sequences, phas24_unique = load_bona_fide_phas(
        args.phas24_loci, args.phas24_products, 24
    )
    unsupported_products, unsupported_intervals, unsupported_sequences, unsupported_summary = (
        build_unsupported_mirna(args, supported_mir_sequences)
    )
    phas_like_products, phas_like_intervals, phas_like_sequences, phas_like_summary = (
        build_phas_like(args, phas21_unique | phas24_unique)
    )
    structural_intervals, structural_sequences, structural_summary = load_structural(args)

    all_intervals = (
        mir_intervals
        + unsupported_intervals
        + phas21_intervals
        + phas24_intervals
        + phas_like_intervals
        + structural_intervals
    )
    genomic_rows = audit_genomic_overlaps(all_intervals)
    genomic_path = args.output_dir / f"{args.release_label}.genomic_interval_overlaps.tsv"
    write_tsv(
        genomic_path,
        genomic_rows,
        [
            "class_a",
            "feature_a",
            "start_a",
            "end_a",
            "strand_a",
            "class_b",
            "feature_b",
            "start_b",
            "end_b",
            "strand_b",
            "chromosome",
            "overlap_start_1based",
            "overlap_end_1based",
            "overlap_nt",
            "strand_relationship",
            "interpretation",
        ],
    )

    short_class_sequences = {
        "miRNA_bona_fide": mir_sequences,
        "miRNA_unsupported_candidate": unsupported_sequences,
        "21-PHAS_bona_fide": phas21_sequences,
        "24-PHAS_bona_fide": phas24_sequences,
        "24-PHAS_like": phas_like_sequences,
    }
    exact_rows = audit_short_exact_conflicts(short_class_sequences)
    exact_path = args.output_dir / f"{args.release_label}.short_product_exact_cross_class_conflicts.tsv"
    write_tsv(
        exact_path,
        exact_rows,
        [
            "mature_sequence_DNA",
            "length_nt",
            "class_count",
            "classes",
            "class_memberships",
            "mapping_interpretation",
        ],
    )
    containment_rows = audit_structural_containment(short_class_sequences, structural_sequences)
    containment_path = args.output_dir / f"{args.release_label}.short_product_structural_containment.tsv"
    write_tsv(
        containment_path,
        containment_rows,
        [
            "mature_sequence_DNA",
            "length_nt",
            "short_RNA_class",
            "short_RNA_member_ids",
            "structural_class",
            "structural_reference_id",
            "orientation_in_structural_reference",
            "mapping_interpretation",
        ],
    )

    interval_counts = Counter(interval.feature_class for interval in all_intervals)
    genomic_pair_counts = Counter(
        "|".join(sorted((str(row["class_a"]), str(row["class_b"])))) for row in genomic_rows
    )
    exact_pair_counts: Counter[str] = Counter()
    for row in exact_rows:
        classes = split_values(str(row["classes"]))
        for left_index, left in enumerate(classes):
            for right in classes[left_index + 1 :]:
                exact_pair_counts[f"{left}|{right}"] += 1
    containment_pair_counts = Counter(
        f"{row['short_RNA_class']}|{row['structural_class']}" for row in containment_rows
    )

    summary_rows: list[dict[str, object]] = []
    for key, value in sorted(unsupported_summary.items()):
        summary_rows.append({"section": "miRNA_unsupported_candidate", "metric": key, "value": value})
    for key, value in sorted(phas_like_summary.items()):
        summary_rows.append({"section": "24_PHAS_like", "metric": key, "value": value})
    for key, value in sorted(structural_summary.items()):
        summary_rows.append({"section": "structural_reference", "metric": key, "value": value})
    for key, value in sorted(interval_counts.items()):
        summary_rows.append({"section": "genomic_features", "metric": key, "value": value})
    for key, value in sorted(genomic_pair_counts.items()):
        summary_rows.append({"section": "genomic_overlap_pairs", "metric": key, "value": value})
    for key, value in sorted(exact_pair_counts.items()):
        summary_rows.append({"section": "exact_short_sequence_conflict_pairs", "metric": key, "value": value})
    for key, value in sorted(containment_pair_counts.items()):
        summary_rows.append({"section": "structural_containment_pairs", "metric": key, "value": value})
    summary_rows.extend(
        [
            {"section": "audit", "metric": "genomic_overlap_rows", "value": len(genomic_rows)},
            {"section": "audit", "metric": "exact_short_sequence_conflict_rows", "value": len(exact_rows)},
            {"section": "audit", "metric": "structural_containment_rows", "value": len(containment_rows)},
            {
                "section": "policy",
                "metric": "MIR_PHAS_coordinate_overlap",
                "value": "review_not_automatic_exclusion",
            },
            {
                "section": "policy",
                "metric": "PHAS_like_scope",
                "value": "113_manually_delimited_24nt_pooled_only_loci_as_exploratory_class",
            },
            {
                "section": "policy",
                "metric": "raw_algorithmic_PHAS_like_scope",
                "value": "excluded_pending_manual_review",
            },
        ]
    )
    summary_path = args.output_dir / f"{args.release_label}.summary.tsv"
    write_tsv(summary_path, summary_rows, ["section", "metric", "value"])

    summary_json_path = args.output_dir / f"{args.release_label}.summary.json"
    with summary_json_path.open("wt", encoding="utf-8") as handle:
        json.dump(
            {
                "release_label": args.release_label,
                "builder_version": VERSION,
                "unsupported_miRNA": unsupported_summary,
                "PHAS_like": phas_like_summary,
                "structural": structural_summary,
                "genomic_overlap_rows": len(genomic_rows),
                "exact_short_sequence_conflict_rows": len(exact_rows),
                "structural_containment_rows": len(containment_rows),
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")

    input_paths = [
        args.mirna_candidate_ledger,
        args.mirna_candidate_arms,
        args.mirna_coordinate_proposals,
        args.mirna_local_alignments,
        args.mirna_bona_fide_loci,
        args.mirna_bona_fide_products,
        args.phas_manual_decisions,
        args.phas21_loci,
        args.phas24_loci,
        args.phas21_products,
        args.phas24_products,
        args.structural_membership,
        args.trna_cca_membership,
    ]
    output_paths = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    manifest_rows: list[dict[str, object]] = []
    command = shlex.join([sys.executable, *sys.argv])
    for path in input_paths:
        manifest_rows.append(
            {
                "record_type": "input",
                "path": str(path.resolve()),
                "sha256": sha256(path),
                "status": "consumed_read_only",
                "assembly_identity": "Zm-B73-REFERENCE-NAM-5.0; B73 RefGen_v5",
                "annotation_identity": "class-specific curated source; see path",
                "command": command,
                "tool_version": VERSION,
            }
        )
    for path in output_paths:
        manifest_rows.append(
            {
                "record_type": "output",
                "path": str(path.resolve()),
                "sha256": sha256(path),
                "status": "current_exploratory_class_or_audit_output",
                "assembly_identity": "Zm-B73-REFERENCE-NAM-5.0; B73 RefGen_v5",
                "annotation_identity": args.release_label,
                "command": command,
                "tool_version": VERSION,
            }
        )
    write_tsv(
        args.manifest_output,
        manifest_rows,
        [
            "record_type",
            "path",
            "sha256",
            "status",
            "assembly_identity",
            "annotation_identity",
            "command",
            "tool_version",
        ],
    )

    print(f"release_label\t{args.release_label}")
    print(f"output_dir\t{args.output_dir}")
    print(f"unsupported_candidate_unique_mapping_sequences\t{unsupported_summary['mapping_safe_unique_sequences']}")
    print(f"manual_24_PHAS_like_loci\t{phas_like_summary['manual_24_PHAS_like_loci']}")
    print(f"genomic_overlap_rows\t{len(genomic_rows)}")
    print(f"exact_short_sequence_conflict_rows\t{len(exact_rows)}")
    print(f"structural_containment_rows\t{len(containment_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
