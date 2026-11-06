#!/usr/bin/env python3
"""Build a reviewable B73 structural-ncRNA mapping reference from completed evidence.

This program makes a preliminary, auditable reference—not a final "gold set".
It keeps every feature-level decision and every collapsed-sequence membership so
that the biological choices can be reviewed before the FASTA files are adopted
by a small-RNA analysis.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from statistics import median
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


SOURCE_CLASSES = ("tRNA", "snRNA", "snoRNA", "rRNA")
OUTPUT_CLASSES = ("tRNA", "snRNA", "snoRNA", "rRNA")
SPLICEOSOMAL_MODELS = {"u1", "u2", "u4", "u5", "u6", "u11", "u12", "u4atac", "u6atac", "u7"}
U3_MODELS = {"u3", "plant_u3", "acea_u3"}
BACTERIAL_OR_ARCHAEAL = re.compile(r"\b(?:bacterial|archaeal)\b", re.IGNORECASE)
UNAMBIGUOUS_DNA = re.compile(r"^[ACGT]+$")
RFAM_RECORD = re.compile(
    r"(?P<target>\S+)\s+RF\d{5}\s+(?P<query>\S+)\s+-\s+"
    r"(?:cm|hmm)\s+\d+\s+\d+\s+\d+\s+\d+\s+[+-]\s+\S+\s+\d+\s+"
    r"\d+(?:\.\d+)?\s+-?\d+(?:\.\d+)?\s+(?P<score>-?\d+(?:\.\d+)?)\s+\S+\s+!"
)


@dataclass
class Candidate:
    candidate_id: str
    source_class: str
    rnacentral_id: str
    description: str
    sequence: str
    source_filter_status: str
    source_filter_reason: str
    maps: list[str] = field(default_factory=list)
    rfam_hits: list[tuple[str, float, str]] = field(default_factory=list)
    selected_class: str = ""
    decision: str = ""
    decision_reason: str = ""

    @property
    def has_exact_b73_mapping(self) -> bool:
        return bool(self.maps)

    @property
    def rfam_models(self) -> list[str]:
        return sorted({hit[0] for hit in self.rfam_hits})

    @property
    def primary_rfam_model(self) -> str:
        if not self.rfam_hits:
            return ""
        return sorted(self.rfam_hits, key=lambda hit: (-hit[1], hit[0]))[0][0]


@dataclass
class TRNACall:
    call_id: str
    sequence: str
    feature_type: str
    locus: str
    isotype: str
    anticodon: str
    score: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = (header, "".join(sequence).upper())
                header, sequence = line[1:], []
            else:
                sequence.append(re.sub(r"\s+", "", line))
    if header is not None:
        records[header.split()[0]] = (header, "".join(sequence).upper())
    return records


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def parse_trnascan_calls(gff_path: Path, fasta_path: Path) -> list[TRNACall]:
    fasta = read_fasta(fasta_path)
    calls: list[TRNACall] = []
    with gff_path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"tRNA", "pseudogene"}:
                continue
            attrs = parse_attributes(fields[8])
            call_id = attrs.get("ID")
            if not call_id or call_id not in fasta:
                raise ValueError(f"{gff_path}: tRNAscan GFF feature lacks matching FASTA sequence: {call_id!r}")
            name = attrs.get("Name", "")
            isotype = attrs.get("isotype", "")
            anticodon = attrs.get("anticodon", "")
            if not isotype or not anticodon:
                match = re.search(r"-([A-Za-z/]+)([ACGT]{3})$", name)
                if match:
                    isotype, anticodon = match.groups()
            locus = f"{fields[0]}:{fields[3]}-{fields[4]}:{fields[6]}"
            calls.append(
                TRNACall(
                    call_id=call_id,
                    sequence=fasta[call_id][1],
                    feature_type=fields[2],
                    locus=locus,
                    isotype=isotype,
                    anticodon=anticodon,
                    score=fields[5],
                )
            )
    if not calls:
        raise ValueError(f"No tRNA/pseudogene features found in {gff_path}")
    return calls


def parse_candidates(membership_path: Path) -> dict[str, Candidate]:
    candidates: dict[str, Candidate] = {}
    with membership_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            candidate = Candidate(
                candidate_id=row["candidate_id"],
                source_class=row["source_class"],
                rnacentral_id=row["rnacentral_id"],
                description=row["source_description"],
                sequence=row["canonical_dna_sequence"],
                source_filter_status=row["source_filter_status"],
                source_filter_reason=row["source_filter_reason"],
            )
            candidates[candidate.candidate_id] = candidate
    return candidates


def add_exact_mappings(candidates: dict[str, Candidate], mapping_dir: Path) -> None:
    for srna_class in SOURCE_CLASSES:
        path = mapping_dir / f"rnacentral_{srna_class}_exact_full_length_B73.sam.gz"
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for raw in handle:
                if raw.startswith("@"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 5 or fields[2] == "*":
                    continue
                candidate = candidates.get(fields[0])
                if candidate is None:
                    raise ValueError(f"{path}: mapping for unknown candidate {fields[0]}")
                strand = "-" if int(fields[1]) & 16 else "+"
                candidate.maps.append(f"{fields[2]}:{fields[3]}:{strand}")
    for candidate in candidates.values():
        candidate.maps = sorted(set(candidate.maps))


def add_rfam_hits(candidates: dict[str, Candidate], table_path: Path) -> list[dict[str, str]]:
    """Read valid tblout records and retain a transparent audit of malformed rows.

    Infernal 1.1.4 produced two concatenated tblout rows in this run. A valid
    record is recognised by its complete ``target RFxxxxx query - cm/hmm ... !``
    structure; an incomplete fragment is not used as evidence.
    """
    anomalies: list[dict[str, str]] = []
    with table_path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("#") or not raw.strip():
                continue
            matches = list(RFAM_RECORD.finditer(raw))
            rf_accessions = len(re.findall(r"RF\d{5}", raw))
            if not matches:
                raise ValueError(f"{table_path}: no complete Rfam record on line {line_number}: {raw!r}")
            if len(matches) != rf_accessions:
                anomalies.append(
                    {
                        "line_number": str(line_number),
                        "RF_accessions_in_raw_line": str(rf_accessions),
                        "complete_records_parsed": str(len(matches)),
                        "decision": "incomplete_concatenated_fragment_excluded; complete_record_retained",
                        "raw_line": raw.rstrip("\n"),
                    }
                )
            for match in matches:
                target, query = match.group("target"), match.group("query")
                if query not in candidates:
                    raise ValueError(f"{table_path}: Rfam hit for unknown candidate {query}")
                candidates[query].rfam_hits.append((target, float(match.group("score")), ""))
    return anomalies


def sno_model(model: str) -> bool:
    name = model.lower()
    return "sno" in name or name in U3_MODELS or name == "u54"


def rrna_model(model: str) -> bool:
    return "rrna" in model.lower()


def assign_candidate_decisions(candidates: dict[str, Candidate]) -> None:
    """Apply the published, conservative review rules to source candidates."""
    for candidate in candidates.values():
        if candidate.source_filter_status != "retain_for_classification":
            candidate.decision = "exclude"
            candidate.decision_reason = candidate.source_filter_reason
            continue
        if candidate.source_class == "tRNA":
            # tRNA sequence selection is governed by B73 tRNAscan-SE locus calls.
            candidate.decision = "defer_to_B73_tRNAscan"
            candidate.decision_reason = "tRNA_reference_uses_nonpseudogene_B73_tRNAscan_sequences"
            continue
        if not candidate.has_exact_b73_mapping:
            candidate.decision = "exclude"
            candidate.decision_reason = "no_exact_full_length_B73_mapping"
            continue
        if not candidate.rfam_hits:
            candidate.decision = "exclude"
            candidate.decision_reason = "no_Rfam_gathering_threshold_hit"
            continue

        models = {model.lower() for model in candidate.rfam_models}
        if candidate.source_class == "snRNA":
            if models & U3_MODELS:
                candidate.selected_class = "snoRNA"
                candidate.decision = "include_reclassified"
                candidate.decision_reason = "Rfam_U3_assignment; retain_as_snoRNA"
            elif models & SPLICEOSOMAL_MODELS:
                candidate.selected_class = "snRNA"
                candidate.decision = "include"
                candidate.decision_reason = "RNAcentral_snRNA_and_spliceosomal_Rfam_assignment_and_exact_B73_mapping"
            else:
                candidate.decision = "exclude"
                candidate.decision_reason = "no_spliceosomal_Rfam_assignment"
        elif candidate.source_class == "snoRNA":
            if any(sno_model(model) for model in models):
                candidate.selected_class = "snoRNA"
                candidate.decision = "include"
                candidate.decision_reason = "RNAcentral_snoRNA_and_snoRNA_Rfam_assignment_and_exact_B73_mapping"
            else:
                candidate.decision = "exclude"
                candidate.decision_reason = "no_snoRNA_Rfam_assignment"
        elif candidate.source_class == "rRNA":
            if BACTERIAL_OR_ARCHAEAL.search(candidate.description):
                candidate.decision = "exclude"
                candidate.decision_reason = "source_description_bacterial_or_archaeal; exclude_from_nuclear_eukaryotic_reference"
            elif any(rrna_model(model) for model in models):
                candidate.selected_class = "rRNA"
                candidate.decision = "include"
                candidate.decision_reason = "RNAcentral_rRNA_and_Rfam_rRNA_assignment_and_exact_B73_mapping"
            else:
                candidate.decision = "exclude"
                candidate.decision_reason = "no_rRNA_Rfam_assignment"


def rrna_family(description: str, primary_model: str) -> str:
    """Assign a reporting family, retaining a traceable Rfam fallback.

    RNAcentral descriptions normally supply the 5S/5.8S/18S/26S family.  Two
    selected B73 records only say "ribosomal RNA"/"partial ribosomal RNA";
    their Rfam models nevertheless distinguish LSU from SSU.  The fallback is
    deliberately limited to these two broad, well-established model groups.
    """
    text = description.lower()
    if re.search(r"\b5s\b", text):
        return "5S"
    if re.search(r"\b5[.]8s\b", text):
        return "5.8S"
    if re.search(r"\b45s\b", text):
        return "45S_precursor"
    if re.search(r"\b(?:18s|ssu|small subunit)\b", text):
        return "18S_SSU"
    if re.search(r"\b(?:26s|28s|lsu|large subunit)\b", text):
        return "26S_LSU"
    model = primary_model.lower()
    if "ssu_rrna" in model:
        return "18S_SSU"
    if "lsu_rrna" in model:
        return "26S_LSU"
    return "unresolved_rRNA"


def create_reference(
    candidates: dict[str, Candidate], trna_calls: list[TRNACall]
) -> tuple[dict[str, dict[str, list[dict[str, str]]]], list[dict[str, str]]]:
    """Return unique sequence memberships by final class and per-feature decisions."""
    references: dict[str, dict[str, list[dict[str, str]]]] = {name: defaultdict(list) for name in OUTPUT_CLASSES}
    decision_rows: list[dict[str, str]] = []

    nonpseudo_trna_sequences = {call.sequence for call in trna_calls if call.feature_type == "tRNA"}
    for call in trna_calls:
        include = call.feature_type == "tRNA"
        decision_rows.append(
            {
                "feature_origin": "tRNAscan-SE",
                "feature_id": call.call_id,
                "source_class": "tRNA",
                "effective_class": "tRNA" if include else "",
                "decision": "include" if include else "exclude",
                "decision_reason": "B73_tRNAscan_nonpseudogene" if include else "tRNAscan_predicted_pseudogene",
                "sequence_length": str(len(call.sequence)),
                "rfam_models": "",
                "exact_B73_loci": call.locus,
                "source_description": f"isotype={call.isotype};anticodon={call.anticodon};score={call.score}",
            }
        )
        if include:
            references["tRNA"][call.sequence].append(
                {
                    "member_origin": "tRNAscan-SE",
                    "member_id": call.call_id,
                    "source_class": "tRNA",
                    "source_description": f"isotype={call.isotype};anticodon={call.anticodon};score={call.score}",
                    "rfam_models": "",
                    "exact_B73_loci": call.locus,
                    "family": f"{call.isotype}-{call.anticodon}" if call.isotype and call.anticodon else "unknown",
                }
            )

    for candidate in candidates.values():
        if candidate.source_class == "tRNA":
            represented = candidate.sequence in nonpseudo_trna_sequences
            candidate.decision = "represented" if represented else "exclude"
            candidate.decision_reason = (
                "sequence_represented_by_B73_tRNAscan_nonpseudogene_reference"
                if represented
                else "not_represented_by_B73_tRNAscan_nonpseudogene_reference"
            )
        decision_rows.append(
            {
                "feature_origin": "RNAcentral",
                "feature_id": candidate.candidate_id,
                "source_class": candidate.source_class,
                "effective_class": candidate.selected_class,
                "decision": candidate.decision,
                "decision_reason": candidate.decision_reason,
                "sequence_length": str(len(candidate.sequence)),
                "rfam_models": ";".join(candidate.rfam_models),
                "exact_B73_loci": ";".join(candidate.maps),
                "source_description": candidate.description,
            }
        )
        if candidate.source_class == "tRNA" and candidate.decision == "represented":
            references["tRNA"][candidate.sequence].append(
                {
                    "member_origin": "RNAcentral",
                    "member_id": candidate.candidate_id,
                    "source_class": candidate.source_class,
                    "source_description": candidate.description,
                    "rfam_models": ";".join(candidate.rfam_models),
                    "exact_B73_loci": ";".join(candidate.maps),
                    "family": candidate.primary_rfam_model or "tRNA",
                }
            )
        elif candidate.decision.startswith("include"):
            family = candidate.primary_rfam_model
            if candidate.selected_class == "rRNA":
                family = rrna_family(candidate.description, candidate.primary_rfam_model)
            references[candidate.selected_class][candidate.sequence].append(
                {
                    "member_origin": "RNAcentral",
                    "member_id": candidate.candidate_id,
                    "source_class": candidate.source_class,
                    "source_description": candidate.description,
                    "rfam_models": ";".join(candidate.rfam_models),
                    "exact_B73_loci": ";".join(candidate.maps),
                    "family": family or "unresolved",
                }
            )
    return references, decision_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curation-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tracked-manifest", type=Path, required=True)
    parser.add_argument(
        "--status",
        choices=("review_required_before_pipeline_use", "approved_for_pipeline_use"),
        default="review_required_before_pipeline_use",
        help="Release status recorded in the manifest.",
    )
    parser.add_argument(
        "--fasta-label",
        default="review",
        help="Version label used in output FASTA filenames (letters, digits, dots, dashes, and underscores only).",
    )
    parser.add_argument(
        "--release-version",
        default="",
        help="Human-readable release version recorded in the manifest, e.g. B73v5-nuclear-structural-ncRNA-v1.",
    )
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", args.fasta_label):
        raise SystemExit("--fasta-label may contain only letters, digits, dots, dashes, and underscores")
    root = args.curation_root
    outdir = args.outdir
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    reference_dir = outdir / "reference"
    reports_dir = outdir / "reports"
    reference_dir.mkdir()
    reports_dir.mkdir()

    membership_path = root / "reports/source_candidate_membership.tsv"
    mapping_dir = root / "genomic_support"
    rfam_table = root / "classification/rfam_candidates.tbl"
    trna_dir = root / "trnascan_genome_evidence"
    trna_gff = trna_dir / "b73v5_trnascan.gff3"
    trna_fasta = trna_dir / "b73v5_trnascan.fa"
    required = (membership_path, rfam_table, trna_gff, trna_fasta)
    absent = [str(path) for path in required if not path.is_file()]
    if absent:
        raise SystemExit("Missing required evidence file(s): " + "; ".join(absent))

    candidates = parse_candidates(membership_path)
    add_exact_mappings(candidates, mapping_dir)
    rfam_anomalies = add_rfam_hits(candidates, rfam_table)
    assign_candidate_decisions(candidates)
    trna_calls = parse_trnascan_calls(trna_gff, trna_fasta)
    references, decision_rows = create_reference(candidates, trna_calls)

    reference_membership_rows: list[dict[str, str]] = []
    reference_family_rows: list[dict[str, str]] = []
    sequence_classes: dict[str, set[str]] = defaultdict(set)
    summary_rows: list[dict[str, str]] = []
    sequence_qc_rows: list[dict[str, str]] = []
    for output_class in OUTPUT_CLASSES:
        sequence_members = references[output_class]
        ordered_sequences = sorted(sequence_members)
        invalid_sequences = [sequence for sequence in ordered_sequences if not UNAMBIGUOUS_DNA.fullmatch(sequence)]
        if invalid_sequences:
            raise SystemExit(
                f"Refusing to write {output_class}: {len(invalid_sequences)} selected sequence(s) contain bases other than A/C/G/T"
            )
        fasta_records: list[tuple[str, str]] = []
        for index, sequence in enumerate(ordered_sequences, start=1):
            reference_id = f"zma_B73v5_{output_class}_{index:05d}"
            fasta_records.append((reference_id, sequence))
            sequence_classes[sequence].add(output_class)
            for member in sequence_members[sequence]:
                reference_membership_rows.append(
                    {
                        "reference_id": reference_id,
                        "class": output_class,
                        "canonical_dna_sequence": sequence,
                        "sequence_length": str(len(sequence)),
                        **member,
                    }
                )
                reference_family_rows.append(
                    {
                        "class": output_class,
                        "reference_id": reference_id,
                        "family": member["family"],
                    }
                )
        write_fasta(reference_dir / f"zma_{output_class}.B73v5.{args.fasta_label}.fa", fasta_records)
        base_counts = Counter("".join(ordered_sequences))
        lengths = [len(sequence) for sequence in ordered_sequences]
        sequence_qc_rows.append(
            {
                "class": output_class,
                "records": str(len(ordered_sequences)),
                "total_bases": str(sum(lengths)),
                "minimum_length": str(min(lengths) if lengths else 0),
                "median_length": str(int(median(lengths)) if lengths else 0),
                "maximum_length": str(max(lengths) if lengths else 0),
                "A_bases": str(base_counts["A"]),
                "C_bases": str(base_counts["C"]),
                "G_bases": str(base_counts["G"]),
                "T_bases": str(base_counts["T"]),
                "N_bases": str(base_counts["N"]),
                "non_ACGT_records": "0",
            }
        )
        summary_rows.append(
            {
                "class": output_class,
                "unique_reference_sequences": str(len(ordered_sequences)),
                "feature_memberships": str(sum(len(members) for members in sequence_members.values())),
                "selection_rule": {
                    "tRNA": "B73 tRNAscan-SE feature type tRNA; predicted pseudogenes excluded; exact duplicate sequences collapsed",
                    "snRNA": "RNAcentral snRNA + exact full-length B73 mapping + spliceosomal Rfam model; exact duplicates collapsed",
                    "snoRNA": "RNAcentral snoRNA + exact full-length B73 mapping + snoRNA Rfam model; includes U3 reclassified from source snRNA; exact duplicates collapsed",
                    "rRNA": "RNAcentral rRNA + exact full-length B73 mapping + rRNA Rfam model; bacterial/archaeal descriptions excluded; exact duplicates collapsed",
                }[output_class],
            }
        )

    write_tsv(
        reports_dir / "feature_decisions.tsv",
        [
            "feature_origin", "feature_id", "source_class", "effective_class", "decision", "decision_reason",
            "sequence_length", "rfam_models", "exact_B73_loci", "source_description",
        ],
        sorted(decision_rows, key=lambda row: (row["feature_origin"], row["source_class"], row["feature_id"])),
    )
    write_tsv(
        reports_dir / "reference_sequence_membership.tsv",
        [
            "reference_id", "class", "canonical_dna_sequence", "sequence_length", "member_origin", "member_id",
            "source_class", "source_description", "rfam_models", "exact_B73_loci", "family",
        ],
        sorted(reference_membership_rows, key=lambda row: (row["class"], row["reference_id"], row["member_origin"], row["member_id"])),
    )
    write_tsv(
        reports_dir / "rfam_table_parse_anomalies.tsv",
        ["line_number", "RF_accessions_in_raw_line", "complete_records_parsed", "decision", "raw_line"],
        rfam_anomalies,
    )
    family_counts = Counter((row["class"], row["family"]) for row in reference_family_rows)
    write_tsv(
        reports_dir / "reference_family_counts.tsv",
        ["class", "family", "feature_memberships"],
        [
            {"class": srna_class, "family": family, "feature_memberships": str(count)}
            for (srna_class, family), count in sorted(family_counts.items())
        ],
    )
    cross_class_rows = []
    for sequence, classes in sorted(sequence_classes.items()):
        if len(classes) > 1:
            cross_class_rows.append(
                {
                    "canonical_dna_sequence": sequence,
                    "sequence_length": str(len(sequence)),
                    "classes": ";".join(sorted(classes)),
                    "decision": "review_required_cross_class_identical_sequence",
                }
            )
    write_tsv(
        reports_dir / "cross_class_exact_sequence_conflicts.tsv",
        ["canonical_dna_sequence", "sequence_length", "classes", "decision"],
        cross_class_rows,
    )
    write_tsv(reports_dir / "reference_class_summary.tsv", list(summary_rows[0]), summary_rows)
    write_tsv(
        reports_dir / "reference_sequence_qc.tsv",
        [
            "class", "records", "total_bases", "minimum_length", "median_length", "maximum_length",
            "A_bases", "C_bases", "G_bases", "T_bases", "N_bases", "non_ACGT_records",
        ],
        sequence_qc_rows,
    )

    manifest_rows = [
        {
            "build_label": outdir.name,
            "release_version": args.release_version,
            "status": args.status,
            "assembly_identity": "Zm-B73-REFERENCE-NAM-5.0; Ensembl Plants release 63",
            "annotation_identity": "RNAcentral v27 typed browser exports; Rfam current covariance models; tRNAscan-SE 2.0.13 B73 nuclear scan",
            "reference_scope": "B73 nuclear structural ncRNA mapping reference; no sample small-RNA sequencing reads used in selection",
            "sequence_representation": "genomic DNA representation (U normalized to T); no post-transcriptional CCA tail added",
            "ambiguous_base_policy": "any candidate with non-ACGT bases is excluded upstream with a decision record; builder fails if a selected FASTA sequence contains a non-ACGT base; no bases are substituted",
            "tRNA_policy": "tRNAscan-SE feature type tRNA retained; predicted pseudogenes excluded; exact duplicate sequences collapsed with locus membership retained",
            "snRNA_policy": "source snRNA + exact B73 mapping + spliceosomal Rfam assignment; U3 reassigned to snoRNA",
            "snoRNA_policy": "source snoRNA + exact B73 mapping + snoRNA Rfam assignment",
            "rRNA_policy": "source rRNA + exact B73 mapping + rRNA Rfam assignment; bacterial/archaeal descriptions excluded",
            "rfam_tblout_parser_policy": "complete target/RF/query/cm-or-hmm records retained; incomplete fragments in concatenated rows excluded and reported",
            "rfam_tblout_anomalous_lines": str(len(rfam_anomalies)),
            "input_membership_path": str(membership_path),
            "input_membership_sha256": sha256(membership_path),
            "input_rfam_table_path": str(rfam_table),
            "input_rfam_table_sha256": sha256(rfam_table),
            "input_trnascan_gff_path": str(trna_gff),
            "input_trnascan_gff_sha256": sha256(trna_gff),
            "input_trnascan_fasta_path": str(trna_fasta),
            "input_trnascan_fasta_sha256": sha256(trna_fasta),
            "builder_script": str(Path(__file__).resolve()),
            "builder_script_sha256": sha256(Path(__file__).resolve()),
        }
    ]
    write_tsv(args.tracked_manifest, list(manifest_rows[0]), manifest_rows)
    status_text = (
        "Approved B73 v5 nuclear structural-ncRNA reference for pipeline use.\n"
        if args.status == "approved_for_pipeline_use"
        else "Preliminary B73 structural-ncRNA reference for review.\nDo not call this a gold set or use it in the pipeline until the report tables are reviewed.\n"
    )
    (outdir / "README.txt").write_text(
        status_text
        + "The reference is independent of this project's sRNA libraries; it derives from B73 v5, RNAcentral, Rfam, and tRNAscan-SE evidence.\n"
        + "Reference FASTA files use genomic DNA representation and deliberately do not add a CCA tail to tRNAs.\n",
        encoding="utf-8",
    )
    print(f"Wrote review reference to {outdir}")
    print(f"Wrote tracked build manifest to {args.tracked_manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
