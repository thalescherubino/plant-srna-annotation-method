#!/usr/bin/env python3
"""Freeze conflict-adjudicated B73-v5 21-PHAS and 24-PHAS references.

The freeze applies explicit user decisions:

* bona-fide MIR loci take precedence over overlapping PHAS calls;
* 24-PHAS loci with phased products exactly contained in curated snRNA or
  snoRNA references are excluded;
* cross-phase locus overlaps retain the phase with the larger maximum
  ``Peak_Howell_score`` across the combined call evidence;
* raw 21-PHAS-like calls are retained as a separate exploratory class without
  manual curation, after removing loci already represented by the adjudicated
  source PHAS annotation; and
* all exclusions and overlap comparisons remain in machine-readable ledgers.

The script filters existing results only. It does not run Phasis or ShortStack.
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
from pathlib import Path
from typing import Iterable


VERSION = "1.0.0"
DNA_RE = re.compile(r"^[ACGT]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--21-reference", dest="reference_21", type=Path, required=True)
    parser.add_argument("--24-reference", dest="reference_24", type=Path, required=True)
    parser.add_argument("--21-products", dest="products_21", type=Path, required=True)
    parser.add_argument("--24-products", dest="products_24", type=Path, required=True)
    parser.add_argument("--21-calls", dest="calls_21", type=Path, required=True)
    parser.add_argument("--24-calls", dest="calls_24", type=Path, required=True)
    parser.add_argument("--cross-class-genomic-overlaps", type=Path, required=True)
    parser.add_argument("--structural-containment", type=Path, required=True)
    parser.add_argument("--21-phas-like-calls", type=Path, action="append", required=True)
    parser.add_argument("--21-phas-like-products", type=Path, action="append", required=True)
    parser.add_argument("--24-phas-like-loci", type=Path, required=True)
    parser.add_argument("--24-phas-like-products", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    parser.add_argument("--release-label", default="Zm_B73v5_PHAS_frozen_v1_2026-08-26")
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


def write_fasta(path: Path, records: Iterable[tuple[str, str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("wt", encoding="utf-8") as handle:
        for identifier, sequence, description in records:
            handle.write(f">{identifier} {description}\n")
            handle.write(sequence + "\n")
            count += 1
    return count


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dna(sequence: str) -> str:
    return sequence.upper().replace("U", "T")


def contig_key(contig: str) -> tuple[int, int | str]:
    match = re.fullmatch(r"chr(\d+)", contig)
    return (0, int(match.group(1))) if match else (1, contig)


def interval_overlap(left: dict[str, object], right: dict[str, object]) -> bool:
    return (
        left["chromosome"] == right["chromosome"]
        and int(left["start_1based"]) <= int(right["end_1based"])
        and int(right["start_1based"]) <= int(left["end_1based"])
    )


def score_evidence(path: Path) -> dict[str, dict[str, object]]:
    values: dict[str, list[float]] = defaultdict(list)
    for row in read_tsv(path):
        if row["Peak_Howell_score"]:
            values[row["identifier"]].append(float(row["Peak_Howell_score"]))
    return {
        identifier: {
            "peak_Howell_score_max": max(scores),
            "Howell_support_call_count": len(scores),
        }
        for identifier, scores in values.items()
    }


def add_exclusion(
    exclusions: dict[tuple[int, str], dict[str, set[str]]],
    phase: int,
    reference_id: str,
    reason: str,
    conflict: str,
) -> None:
    record = exclusions[(phase, reference_id)]
    record["reasons"].add(reason)
    if conflict:
        record["conflicts"].add(conflict)


def parse_conflicts(
    args: argparse.Namespace,
    refs: dict[int, dict[str, dict[str, str]]],
    scores: dict[int, dict[str, dict[str, object]]],
) -> tuple[dict[tuple[int, str], dict[str, set[str]]], list[dict[str, object]]]:
    exclusions: dict[tuple[int, str], dict[str, set[str]]] = defaultdict(
        lambda: {"reasons": set(), "conflicts": set()}
    )
    overlap_adjudications: list[dict[str, object]] = []

    for row in read_tsv(args.cross_class_genomic_overlaps):
        classes = {row["class_a"], row["class_b"]}
        if "miRNA_bona_fide" in classes:
            for phase in (21, 24):
                phas_class = f"{phase}-PHAS_bona_fide"
                if phas_class not in classes:
                    continue
                phas_id = row["feature_a"] if row["class_a"] == phas_class else row["feature_b"]
                mir_id = row["feature_a"] if row["class_a"] == "miRNA_bona_fide" else row["feature_b"]
                add_exclusion(
                    exclusions,
                    phase,
                    phas_id,
                    "overlaps_bona_fide_MIR_locus",
                    mir_id,
                )

        if classes == {"21-PHAS_bona_fide", "24-PHAS_bona_fide"}:
            id21 = row["feature_a"] if row["class_a"] == "21-PHAS_bona_fide" else row["feature_b"]
            id24 = row["feature_a"] if row["class_a"] == "24-PHAS_bona_fide" else row["feature_b"]
            source21 = refs[21][id21]["source_identifier"]
            source24 = refs[24][id24]["source_identifier"]
            score21 = scores[21].get(source21, {}).get("peak_Howell_score_max")
            score24 = scores[24].get(source24, {}).get("peak_Howell_score_max")
            if score21 is None or score24 is None:
                raise ValueError(
                    f"Cross-phase overlap lacks Howell evidence: {id21}={score21}; {id24}={score24}"
                )
            if float(score21) == float(score24):
                raise ValueError(f"Cross-phase Howell-score tie requires explicit adjudication: {id21}; {id24}")
            winner_phase, winner_id, winner_score, loser_phase, loser_id, loser_score = (
                (21, id21, score21, 24, id24, score24)
                if float(score21) > float(score24)
                else (24, id24, score24, 21, id21, score21)
            )
            add_exclusion(
                exclusions,
                loser_phase,
                loser_id,
                "lower_maximum_Peak_Howell_score_in_21_24_overlap",
                winner_id,
            )
            overlap_adjudications.append(
                {
                    "chromosome": row["chromosome"],
                    "overlap_start_1based": row["overlap_start_1based"],
                    "overlap_end_1based": row["overlap_end_1based"],
                    "overlap_nt": row["overlap_nt"],
                    "21_reference_locus_id": id21,
                    "21_source_identifier": source21,
                    "21_maximum_Peak_Howell_score": f"{float(score21):.12g}",
                    "24_reference_locus_id": id24,
                    "24_source_identifier": source24,
                    "24_maximum_Peak_Howell_score": f"{float(score24):.12g}",
                    "retained_phase_nt": winner_phase,
                    "retained_reference_locus_id": winner_id,
                    "retained_maximum_Peak_Howell_score": f"{float(winner_score):.12g}",
                    "excluded_phase_nt": loser_phase,
                    "excluded_reference_locus_id": loser_id,
                    "excluded_maximum_Peak_Howell_score": f"{float(loser_score):.12g}",
                    "decision_rule": "retain_locus_with_higher_maximum_Peak_Howell_score",
                }
            )

    for row in read_tsv(args.structural_containment):
        if row["short_RNA_class"] != "24-PHAS_bona_fide":
            continue
        if row["structural_class"] not in {"snRNA", "snoRNA"}:
            continue
        for locus_id in row["short_RNA_member_ids"].split(","):
            add_exclusion(
                exclusions,
                24,
                locus_id,
                f"phased_product_exactly_contained_in_curated_{row['structural_class']}",
                row["structural_reference_id"],
            )

    unknown = [key for key in exclusions if key[1] not in refs[key[0]]]
    if unknown:
        raise ValueError(f"Exclusion references unknown PHAS loci: {unknown[:5]}")
    return exclusions, overlap_adjudications


def freeze_phase(
    args: argparse.Namespace,
    phase: int,
    source_rows: list[dict[str, str]],
    product_rows: list[dict[str, str]],
    score_rows: dict[str, dict[str, object]],
    exclusions: dict[tuple[int, str], dict[str, set[str]]],
) -> tuple[list[dict[str, object]], list[dict[str, str]], list[dict[str, object]], set[str]]:
    retained_source = [
        row for row in source_rows if (phase, row["reference_locus_id"]) not in exclusions
    ]
    retained_source.sort(
        key=lambda row: (
            contig_key(row["chromosome"]),
            int(row["start_1based"]),
            int(row["end_1based"]),
            row["reference_locus_id"],
        )
    )
    frozen_rows: list[dict[str, object]] = []
    source_to_frozen: dict[str, str] = {}
    reference_to_frozen: dict[str, str] = {}
    for index, row in enumerate(retained_source, start=1):
        frozen_id = f"ZmaB73v5-{phase}PHAS-FRZ1-{index:06d}"
        source_identifier = row["source_identifier"]
        evidence = score_rows.get(source_identifier, {})
        record: dict[str, object] = dict(row)
        record.update(
            {
                "frozen_locus_id": frozen_id,
                "freeze_release": args.release_label,
                "freeze_status": "retained_canonical_PHAS",
                "maximum_Peak_Howell_score": (
                    f"{float(evidence['peak_Howell_score_max']):.12g}"
                    if evidence
                    else ""
                ),
                "Howell_support_call_count": evidence.get("Howell_support_call_count", ""),
            }
        )
        frozen_rows.append(record)
        reference_to_frozen[row["reference_locus_id"]] = frozen_id
        if row["source_method"] == "Phasis":
            source_to_frozen[source_identifier] = frozen_id

    retained_products = [row for row in product_rows if row["identifier"] in source_to_frozen]
    retained_product_loci = {row["identifier"] for row in retained_products}
    if retained_product_loci != set(source_to_frozen):
        missing = sorted(set(source_to_frozen) - retained_product_loci)
        raise ValueError(f"Retained {phase}-PHAS Phasis loci lack products: {missing[:5]}")

    frozen_fields = list(source_rows[0]) + [
        "frozen_locus_id",
        "freeze_release",
        "freeze_status",
        "maximum_Peak_Howell_score",
        "Howell_support_call_count",
    ]
    write_tsv(
        args.output_dir / f"{phase}_PHAS.frozen.loci.tsv",
        frozen_rows,
        frozen_fields,
    )
    write_tsv(
        args.output_dir / f"{phase}_PHAS.frozen.phasiRNAs.tsv",
        retained_products,
        list(product_rows[0]),
    )
    with (args.output_dir / f"{phase}_PHAS.frozen.loci.bed").open("wt", encoding="utf-8") as handle:
        for row in frozen_rows:
            handle.write(
                "\t".join(
                    [
                        str(row["chromosome"]),
                        str(int(row["start_1based"]) - 1),
                        str(row["end_1based"]),
                        str(row["frozen_locus_id"]),
                        str(row["maximum_Peak_Howell_score"] or 0),
                        ".",
                    ]
                )
                + "\n"
            )

    sequence_members: dict[str, set[str]] = defaultdict(set)
    for row in retained_products:
        sequence = dna(row["tag_seq"])
        sequence_members[sequence].add(source_to_frozen[row["identifier"]])
    sequence_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str, str]] = []
    for index, sequence in enumerate(sorted(sequence_members), start=1):
        sequence_id = f"ZmaB73v5-{phase}PHAS-FRZ1-SEQ{index:07d}"
        members = sorted(sequence_members[sequence])
        sequence_rows.append(
            {
                "sequence_reference_id": sequence_id,
                "sequence_DNA": sequence,
                "length_nt": len(sequence),
                "frozen_locus_count": len(members),
                "frozen_locus_ids": ",".join(members),
            }
        )
        fasta_records.append(
            (sequence_id, sequence, f"class={phase}_PHAS frozen_locus_count={len(members)}")
        )
    write_tsv(
        args.output_dir / f"{phase}_PHAS.frozen.sequence_equivalence.tsv",
        sequence_rows,
        [
            "sequence_reference_id",
            "sequence_DNA",
            "length_nt",
            "frozen_locus_count",
            "frozen_locus_ids",
        ],
    )
    write_fasta(
        args.output_dir / f"{phase}_PHAS.frozen.sequence_equivalence.fa",
        fasta_records,
    )

    mappings = [
        {
            "phase_nt": phase,
            "frozen_locus_id": reference_to_frozen[row["reference_locus_id"]],
            "source_reference_locus_id": row["reference_locus_id"],
            "source_identifier": row["source_identifier"],
            "source_method": row["source_method"],
        }
        for row in retained_source
    ]
    return frozen_rows, retained_products, mappings, set(sequence_members)


def parse_identifier_interval(identifier: str) -> dict[str, object]:
    match = re.fullmatch(r"([^:]+):(\d+)\.\.(\d+)", identifier)
    if not match:
        raise ValueError(f"Unrecognized PHAS identifier: {identifier!r}")
    contig, start, end = match.groups()
    return {
        "chromosome": contig,
        "start_1based": int(start),
        "end_1based": int(end),
    }


def build_21_phas_like(
    args: argparse.Namespace,
    initial_21: list[dict[str, str]],
    final_phas_sequences: set[str],
) -> dict[str, int]:
    call_evidence: dict[str, dict[str, object]] = defaultdict(
        lambda: {"modes": set(), "scores": [], "rows": 0}
    )
    for mode, path in zip(("nonpooled", "pooled"), args.__dict__["21_phas_like_calls"]):
        for row in read_tsv(path):
            record = call_evidence[row["identifier"]]
            record["modes"].add(mode)
            record["rows"] = int(record["rows"]) + 1
            if row["Peak_Howell_score"]:
                record["scores"].append(float(row["Peak_Howell_score"]))

    initial_intervals = [
        {
            "reference_locus_id": row["reference_locus_id"],
            "chromosome": row["chromosome"],
            "start_1based": int(row["start_1based"]),
            "end_1based": int(row["end_1based"]),
        }
        for row in initial_21
    ]
    source_ledger: list[dict[str, object]] = []
    eligible: list[dict[str, object]] = []
    for identifier in sorted(call_evidence):
        coords = parse_identifier_interval(identifier)
        overlaps = [row["reference_locus_id"] for row in initial_intervals if interval_overlap(coords, row)]
        evidence = call_evidence[identifier]
        status = (
            "represented_by_source_bona_fide_21PHAS_annotation"
            if overlaps
            else "include_as_distinct_exploratory_21_PHAS_like"
        )
        row = {
            "source_identifier": identifier,
            **coords,
            "evidence_modes": ",".join(sorted(evidence["modes"])),
            "source_call_row_count": evidence["rows"],
            "maximum_Peak_Howell_score": (
                f"{max(evidence['scores']):.12g}" if evidence["scores"] else ""
            ),
            "overlapping_source_bona_fide_21PHAS_loci": ",".join(sorted(overlaps)),
            "class_assignment": status,
            "manual_curation_status": "not_required_by_user_decision",
        }
        source_ledger.append(row)
        if not overlaps:
            eligible.append(row)

    eligible.sort(
        key=lambda row: (
            contig_key(str(row["chromosome"])),
            int(row["start_1based"]),
            int(row["end_1based"]),
        )
    )
    identifier_to_locus: dict[str, str] = {}
    loci: list[dict[str, object]] = []
    for index, row in enumerate(eligible, start=1):
        locus_id = f"ZmaB73v5-21PHASLIKE-FRZ1-{index:06d}"
        identifier_to_locus[str(row["source_identifier"])] = locus_id
        loci.append(
            {
                "phas_like_locus_id": locus_id,
                **row,
                "phase_nt": 21,
                "reference_status": "exploratory_not_bona_fide_PHAS",
                "product_semantics": "Phasis_reported_product_candidate_not_validated_phasiRNA",
            }
        )
    write_tsv(
        args.output_dir / "21_PHAS_like.source_call_ledger.tsv",
        source_ledger,
        [
            "source_identifier",
            "chromosome",
            "start_1based",
            "end_1based",
            "evidence_modes",
            "source_call_row_count",
            "maximum_Peak_Howell_score",
            "overlapping_source_bona_fide_21PHAS_loci",
            "class_assignment",
            "manual_curation_status",
        ],
    )
    write_tsv(
        args.output_dir / "21_PHAS_like.exploratory.loci.tsv",
        loci,
        [
            "phas_like_locus_id",
            "source_identifier",
            "chromosome",
            "start_1based",
            "end_1based",
            "phase_nt",
            "evidence_modes",
            "source_call_row_count",
            "maximum_Peak_Howell_score",
            "overlapping_source_bona_fide_21PHAS_loci",
            "class_assignment",
            "manual_curation_status",
            "reference_status",
            "product_semantics",
        ],
    )

    aggregate: dict[tuple[str, str], dict[str, object]] = {}
    raw_product_rows = 0
    for mode, path in zip(("nonpooled", "pooled"), args.__dict__["21_phas_like_products"]):
        for row in read_tsv(path):
            identifier = row["identifier"]
            if identifier not in identifier_to_locus:
                continue
            raw_product_rows += 1
            sequence = dna(row["tag_seq"])
            if not DNA_RE.fullmatch(sequence) or len(sequence) != 21:
                raise ValueError(f"Invalid 21-PHAS-like product: {sequence!r}")
            locus_id = identifier_to_locus[identifier]
            record = aggregate.setdefault(
                (locus_id, sequence),
                {
                    "phas_like_locus_id": locus_id,
                    "sequence_DNA": sequence,
                    "source_row_count": 0,
                    "abundance_sum": 0.0,
                    "abundance_max": 0.0,
                    "source_modes": set(),
                    "source_libraries": set(),
                    "register_classes": set(),
                    "source_identifiers": set(),
                },
            )
            abundance = float(row["abun"])
            record["source_row_count"] = int(record["source_row_count"]) + 1
            record["abundance_sum"] = float(record["abundance_sum"]) + abundance
            record["abundance_max"] = max(float(record["abundance_max"]), abundance)
            record["source_modes"].add(mode)
            record["source_libraries"].add(row["alib"])
            record["register_classes"].add(row["register_class"])
            record["source_identifiers"].add(identifier)
    product_rows: list[dict[str, object]] = []
    for index, key in enumerate(sorted(aggregate), start=1):
        row = aggregate[key]
        product_rows.append(
            {
                "phas_like_product_id": f"ZmaB73v5-21PHASLIKE-FRZ1-P{index:06d}",
                "phas_like_locus_id": row["phas_like_locus_id"],
                "sequence_DNA": row["sequence_DNA"],
                "length_nt": 21,
                "source_row_count": row["source_row_count"],
                "reported_abundance_sum": f"{float(row['abundance_sum']):.12g}",
                "reported_abundance_max": f"{float(row['abundance_max']):.12g}",
                "source_modes": ",".join(sorted(row["source_modes"])),
                "source_libraries": ",".join(sorted(row["source_libraries"])),
                "register_classes": ",".join(sorted(row["register_classes"])),
                "source_identifiers": ",".join(sorted(row["source_identifiers"])),
                "identical_to_frozen_bona_fide_PHAS_product": (
                    "yes" if row["sequence_DNA"] in final_phas_sequences else "no"
                ),
                "product_status": "21_PHAS_like_product_candidate_not_validated_phasiRNA",
            }
        )
    product_fields = [
        "phas_like_product_id",
        "phas_like_locus_id",
        "sequence_DNA",
        "length_nt",
        "source_row_count",
        "reported_abundance_sum",
        "reported_abundance_max",
        "source_modes",
        "source_libraries",
        "register_classes",
        "source_identifiers",
        "identical_to_frozen_bona_fide_PHAS_product",
        "product_status",
    ]
    write_tsv(args.output_dir / "21_PHAS_like.exploratory.product_candidates.tsv", product_rows, product_fields)

    by_sequence: dict[str, set[str]] = defaultdict(set)
    for row in product_rows:
        if row["identical_to_frozen_bona_fide_PHAS_product"] == "no":
            by_sequence[str(row["sequence_DNA"])].add(str(row["phas_like_locus_id"]))
    sequence_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str, str]] = []
    for index, sequence in enumerate(sorted(by_sequence), start=1):
        sequence_id = f"ZmaB73v5-21PHASLIKE-FRZ1-SEQ{index:06d}"
        loci_for_sequence = sorted(by_sequence[sequence])
        sequence_rows.append(
            {
                "sequence_reference_id": sequence_id,
                "sequence_DNA": sequence,
                "length_nt": 21,
                "locus_count": len(loci_for_sequence),
                "phas_like_locus_ids": ",".join(loci_for_sequence),
            }
        )
        fasta_records.append(
            (sequence_id, sequence, f"class=21_PHAS_like locus_count={len(loci_for_sequence)}")
        )
    write_tsv(
        args.output_dir / "21_PHAS_like.exploratory.sequence_equivalence.tsv",
        sequence_rows,
        ["sequence_reference_id", "sequence_DNA", "length_nt", "locus_count", "phas_like_locus_ids"],
    )
    write_fasta(args.output_dir / "21_PHAS_like.exploratory.sequence_equivalence.fa", fasta_records)
    return {
        "raw_source_identifiers": len(source_ledger),
        "source_identifiers_already_represented_or_overlapping_bona_fide_21PHAS": (
            len(source_ledger) - len(loci)
        ),
        "distinct_exploratory_loci": len(loci),
        "source_product_rows_for_distinct_loci": raw_product_rows,
        "locus_sequence_product_candidates": len(product_rows),
        "mapping_safe_unique_sequences": len(sequence_rows),
    }


def rebuild_24_phas_like(
    args: argparse.Namespace,
    final_phas_sequences: set[str],
) -> dict[str, int]:
    loci = read_tsv(args.__dict__["24_phas_like_loci"])
    products = read_tsv(args.__dict__["24_phas_like_products"])
    by_sequence: dict[str, set[str]] = defaultdict(set)
    conflicting_products = 0
    for row in products:
        sequence = dna(row["mature_sequence_DNA"])
        if sequence in final_phas_sequences:
            conflicting_products += 1
            continue
        by_sequence[sequence].add(row["phas_like_locus_id"])
    sequence_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str, str]] = []
    for index, sequence in enumerate(sorted(by_sequence), start=1):
        sequence_id = f"ZmaB73v5-24PHASLIKE-FRZ1-SEQ{index:07d}"
        locus_ids = sorted(by_sequence[sequence])
        sequence_rows.append(
            {
                "sequence_reference_id": sequence_id,
                "sequence_DNA": sequence,
                "length_nt": 24,
                "locus_count": len(locus_ids),
                "phas_like_locus_ids": ",".join(locus_ids),
            }
        )
        fasta_records.append(
            (sequence_id, sequence, f"class=24_PHAS_like locus_count={len(locus_ids)}")
        )
    write_tsv(
        args.output_dir / "24_PHAS_like.exploratory.sequence_equivalence.tsv",
        sequence_rows,
        ["sequence_reference_id", "sequence_DNA", "length_nt", "locus_count", "phas_like_locus_ids"],
    )
    write_fasta(args.output_dir / "24_PHAS_like.exploratory.sequence_equivalence.fa", fasta_records)
    return {
        "manual_exploratory_loci": len(loci),
        "locus_sequence_product_candidates": len(products),
        "products_identical_to_frozen_bona_fide_PHAS": conflicting_products,
        "mapping_safe_unique_sequences": len(sequence_rows),
    }


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    source_21 = read_tsv(args.reference_21)
    source_24 = read_tsv(args.reference_24)
    products_21 = read_tsv(args.products_21)
    products_24 = read_tsv(args.products_24)
    refs = {
        21: {row["reference_locus_id"]: row for row in source_21},
        24: {row["reference_locus_id"]: row for row in source_24},
    }
    scores = {21: score_evidence(args.calls_21), 24: score_evidence(args.calls_24)}
    exclusions, overlap_adjudications = parse_conflicts(args, refs, scores)

    frozen21, frozen_products21, mappings21, sequences21 = freeze_phase(
        args, 21, source_21, products_21, scores[21], exclusions
    )
    frozen24, frozen_products24, mappings24, sequences24 = freeze_phase(
        args, 24, source_24, products_24, scores[24], exclusions
    )
    final_phas_sequences = sequences21 | sequences24

    write_tsv(
        args.output_dir / "frozen_locus_id_map.tsv",
        [*mappings21, *mappings24],
        ["phase_nt", "frozen_locus_id", "source_reference_locus_id", "source_identifier", "source_method"],
    )
    write_tsv(
        args.output_dir / "21_24_overlap_Howell_adjudication.tsv",
        sorted(
            overlap_adjudications,
            key=lambda row: (contig_key(str(row["chromosome"])), int(row["overlap_start_1based"])),
        ),
        [
            "chromosome",
            "overlap_start_1based",
            "overlap_end_1based",
            "overlap_nt",
            "21_reference_locus_id",
            "21_source_identifier",
            "21_maximum_Peak_Howell_score",
            "24_reference_locus_id",
            "24_source_identifier",
            "24_maximum_Peak_Howell_score",
            "retained_phase_nt",
            "retained_reference_locus_id",
            "retained_maximum_Peak_Howell_score",
            "excluded_phase_nt",
            "excluded_reference_locus_id",
            "excluded_maximum_Peak_Howell_score",
            "decision_rule",
        ],
    )

    exclusion_rows: list[dict[str, object]] = []
    for (phase, reference_id), details in sorted(
        exclusions.items(),
        key=lambda item: (
            item[0][0],
            contig_key(refs[item[0][0]][item[0][1]]["chromosome"]),
            int(refs[item[0][0]][item[0][1]]["start_1based"]),
        ),
    ):
        source = refs[phase][reference_id]
        evidence = scores[phase].get(source["source_identifier"], {})
        exclusion_rows.append(
            {
                "phase_nt": phase,
                "reference_locus_id": reference_id,
                "source_identifier": source["source_identifier"],
                "chromosome": source["chromosome"],
                "start_1based": source["start_1based"],
                "end_1based": source["end_1based"],
                "maximum_Peak_Howell_score": (
                    f"{float(evidence['peak_Howell_score_max']):.12g}" if evidence else ""
                ),
                "exclusion_reasons": ";".join(sorted(details["reasons"])),
                "conflicting_features": ";".join(sorted(details["conflicts"])),
                "frozen_annotation_status": "excluded_from_canonical_PHAS",
                "DA_follow_up": "retain_in_conflict_ledger_for_post_DA_cross_reference",
            }
        )
    write_tsv(
        args.output_dir / "PHAS_exclusion_ledger.tsv",
        exclusion_rows,
        [
            "phase_nt",
            "reference_locus_id",
            "source_identifier",
            "chromosome",
            "start_1based",
            "end_1based",
            "maximum_Peak_Howell_score",
            "exclusion_reasons",
            "conflicting_features",
            "frozen_annotation_status",
            "DA_follow_up",
        ],
    )

    like21_summary = build_21_phas_like(args, source_21, final_phas_sequences)
    like24_summary = rebuild_24_phas_like(args, final_phas_sequences)

    reason_counts: Counter[str] = Counter()
    for details in exclusions.values():
        for reason in details["reasons"]:
            reason_counts[reason] += 1
    summary = {
        "release_label": args.release_label,
        "source_21_PHAS_loci": len(source_21),
        "frozen_21_PHAS_loci": len(frozen21),
        "excluded_21_PHAS_loci": len(source_21) - len(frozen21),
        "source_24_PHAS_loci": len(source_24),
        "frozen_24_PHAS_loci": len(frozen24),
        "excluded_24_PHAS_loci": len(source_24) - len(frozen24),
        "frozen_21_PHAS_product_rows": len(frozen_products21),
        "frozen_24_PHAS_product_rows": len(frozen_products24),
        "frozen_21_PHAS_unique_sequences": len(sequences21),
        "frozen_24_PHAS_unique_sequences": len(sequences24),
        "MIR_overlapping_PHAS_loci_excluded": sum(
            "overlaps_bona_fide_MIR_locus" in details["reasons"] for details in exclusions.values()
        ),
        "structural_product_24_PHAS_loci_excluded": sum(
            any(reason.startswith("phased_product_exactly_contained") for reason in details["reasons"])
            for details in exclusions.values()
        ),
        "cross_phase_overlap_pairs_adjudicated": len(overlap_adjudications),
        "21_PHAS_like": like21_summary,
        "24_PHAS_like": like24_summary,
        "exclusion_reason_counts": dict(sorted(reason_counts.items())),
    }
    with (args.output_dir / "freeze_summary.json").open("wt", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    summary_rows: list[dict[str, object]] = []
    for key, value in summary.items():
        if isinstance(value, dict):
            for nested_key, nested_value in value.items():
                summary_rows.append({"section": key, "metric": nested_key, "value": nested_value})
        else:
            summary_rows.append({"section": "freeze", "metric": key, "value": value})
    write_tsv(args.output_dir / "freeze_summary.tsv", summary_rows, ["section", "metric", "value"])

    input_paths = [
        args.reference_21,
        args.reference_24,
        args.products_21,
        args.products_24,
        args.calls_21,
        args.calls_24,
        args.cross_class_genomic_overlaps,
        args.structural_containment,
        *args.__dict__["21_phas_like_calls"],
        *args.__dict__["21_phas_like_products"],
        args.__dict__["24_phas_like_loci"],
        args.__dict__["24_phas_like_products"],
    ]
    output_paths = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    command = shlex.join([sys.executable, *sys.argv])
    manifest_rows: list[dict[str, object]] = []
    for path in input_paths:
        manifest_rows.append(
            {
                "record_type": "input",
                "path": str(path.resolve()),
                "sha256": sha256(path),
                "assembly_identity": "Zm-B73-REFERENCE-NAM-5.0; B73 RefGen_v5",
                "annotation_identity": "source curated PHAS or conflict audit; see path",
                "command": command,
                "tool_version": VERSION,
                "status": "consumed_read_only",
            }
        )
    for path in output_paths:
        manifest_rows.append(
            {
                "record_type": "output",
                "path": str(path.resolve()),
                "sha256": sha256(path),
                "assembly_identity": "Zm-B73-REFERENCE-NAM-5.0; B73 RefGen_v5",
                "annotation_identity": args.release_label,
                "command": command,
                "tool_version": VERSION,
                "status": "frozen_current_PHAS_reference_or_exploratory_PHAS_like",
            }
        )
    write_tsv(
        args.manifest_output,
        manifest_rows,
        [
            "record_type",
            "path",
            "sha256",
            "assembly_identity",
            "annotation_identity",
            "command",
            "tool_version",
            "status",
        ],
    )

    print(f"release_label\t{args.release_label}")
    print(f"frozen_21_PHAS_loci\t{len(frozen21)}")
    print(f"frozen_24_PHAS_loci\t{len(frozen24)}")
    print(f"excluded_PHAS_loci\t{len(exclusion_rows)}")
    print(f"exploratory_21_PHAS_like_loci\t{like21_summary['distinct_exploratory_loci']}")
    print(f"exploratory_24_PHAS_like_loci\t{like24_summary['manual_exploratory_loci']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
