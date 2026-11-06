#!/usr/bin/env python3
"""Build the strict B73 RefGen_v5 bona-fide MIR locus and product release.

The builder combines three evidence branches:

1. source precursors that exactly match B73 v5 and passed the original miRScore run;
2. coordinate-remediated B73-v5 models that passed their required miRScore rerun;
3. originally input-incompatible models repaired from B73 v5 that passed rerunning.

Only passing primary models or passing ``-rescue`` alternatives enter the release.
Same-strand overlapping precursor models are represented by one genomic MIR locus,
while distinct genomic loci are retained even when they produce identical products.
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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable
from urllib.parse import quote


VERSION = "1.0.0"
ASSEMBLY = "Zm-B73-REFERENCE-NAM-5.0; B73 RefGen_v5"
ANNOTATION = (
    "strict bona-fide MIR release from miRBase v23, miRScore-study Supplementary "
    "File S4, and the expanded B73-v5 miR2275 set; miRScore 0.3.7 validation"
)
DNA_RE = re.compile(r"^[ACGT]+$")


@dataclass
class ProductEvidence:
    role: str
    arm: str
    sequence: str
    relative_start: int
    relative_end: int
    genomic_start: int = 0
    genomic_end: int = 0
    position_resolution: str = "miRScore_reported_coordinates"


@dataclass
class AcceptedModel:
    model_id: str
    branch: str
    source_internal_id: str
    display_name: str
    source: str
    original_id: str
    predictor: str
    family: str
    chromosome: str
    start: int
    end: int
    strand: str
    coordinate_method: str
    precursor: str
    evidence_type: str
    mirscore_version: str
    total_reads: int
    supporting_libraries: int
    precision: float
    products: list[ProductEvidence] = field(default_factory=list)
    release_locus_id: str = ""
    release_transcript_id: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-validation", type=Path, required=True)
    parser.add_argument("--coordinate-proposals", type=Path, required=True)
    parser.add_argument("--remediation-validation", type=Path, required=True)
    parser.add_argument("--remediation-products", type=Path, required=True)
    parser.add_argument("--input-remediation-validation", type=Path, required=True)
    parser.add_argument("--input-remediation-products", type=Path, required=True)
    parser.add_argument("--release-policy", type=Path, required=True)
    parser.add_argument("--release-decisions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    parser.add_argument("--release-prefix", default="Zm_B73v5_bona_fide_miRNA_v1")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
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


def dna(sequence: str) -> str:
    return sequence.upper().replace("U", "T")


def integer(value: str | int | float) -> int:
    if value in ("", None):
        return 0
    return int(float(value))


def floating(value: str | int | float) -> float:
    if value in ("", None):
        return 0.0
    return float(value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def natural_contig_key(contig: str) -> tuple[int, int | str]:
    match = re.fullmatch(r"chr(\d+)", contig, flags=re.IGNORECASE)
    return (0, int(match.group(1))) if match else (1, contig)


def original_evidence_type(row: dict[str, str]) -> str:
    if row.get("mirscore_result") == "Pass":
        return "primary"
    if row.get("mirscore_rescue_result") == "Pass":
        return "rescue"
    return ""


def products_from_validation(row: dict[str, str], evidence_type: str) -> list[ProductEvidence]:
    prefix = "mirscore_" if evidence_type == "primary" else "mirscore_rescue_"
    products: list[ProductEvidence] = []
    for short_role, role in (("m", "mature"), ("ms", "star")):
        sequence = dna(row.get(f"{prefix}{short_role}Seq", ""))
        start = integer(row.get(f"{prefix}{short_role}Start", ""))
        end = integer(row.get(f"{prefix}{short_role}Stop", ""))
        if not sequence or not start or not end:
            raise ValueError(
                f"Passing candidate {row.get('internal_id', '')} lacks {role} product coordinates"
            )
        products.append(ProductEvidence(role, "", sequence, start, end))
    for index, product in enumerate(sorted(products, key=lambda item: item.relative_start)):
        product.arm = "5p" if index == 0 else "3p"
    return products


def add_genomic_product_coordinates(model: AcceptedModel) -> None:
    if model.start < 1 or model.end < model.start:
        raise ValueError(f"Invalid model interval for {model.model_id}")
    if len(model.precursor) != model.end - model.start + 1:
        raise ValueError(
            f"Precursor length/interval mismatch for {model.model_id}: "
            f"{len(model.precursor)} versus {model.end - model.start + 1}"
        )
    if not DNA_RE.fullmatch(model.precursor):
        raise ValueError(f"Non-ACGT precursor in accepted model {model.model_id}")
    if len(model.products) != 2:
        raise ValueError(f"Accepted model {model.model_id} does not have two products")
    for product in model.products:
        if not 1 <= product.relative_start <= product.relative_end <= len(model.precursor):
            raise ValueError(f"Invalid product interval in {model.model_id}")
        observed = model.precursor[product.relative_start - 1 : product.relative_end]
        if observed != product.sequence:
            starts = [
                match.start()
                for match in re.finditer(f"(?={re.escape(product.sequence)})", model.precursor)
            ]
            if len(starts) != 1:
                raise ValueError(
                    f"Product sequence does not uniquely match the B73-v5 precursor for "
                    f"{model.model_id}: sequence={product.sequence}; matches={len(starts)}"
                )
            product.relative_start = starts[0] + 1
            product.relative_end = starts[0] + len(product.sequence)
            product.position_resolution = "unique_precursor_sequence_match_corrected_miRScore_offset"
        if not 19 <= len(product.sequence) <= 24:
            raise ValueError(f"Product length outside 19-24 nt in {model.model_id}")
        if model.strand == "+":
            product.genomic_start = model.start + product.relative_start - 1
            product.genomic_end = model.start + product.relative_end - 1
        elif model.strand == "-":
            product.genomic_start = model.end - product.relative_end + 1
            product.genomic_end = model.end - product.relative_start + 1
        else:
            raise ValueError(f"Invalid strand in {model.model_id}: {model.strand}")


def exact_models(
    candidates: dict[str, dict[str, str]], proposals: list[dict[str, str]]
) -> list[AcceptedModel]:
    models: list[AcceptedModel] = []
    for proposal in proposals:
        if proposal.get("coordinate_acceptance") != "accepted":
            continue
        if proposal.get("coordinate_method") != "bowtie_exact":
            continue
        candidate = candidates[proposal["internal_id"]]
        evidence = original_evidence_type(candidate)
        if not evidence:
            continue
        prefix = "mirscore_" if evidence == "primary" else "mirscore_rescue_"
        model = AcceptedModel(
            model_id=f"EXACT_{proposal['projection_id']}",
            branch="original_exact_B73v5",
            source_internal_id=candidate["internal_id"],
            display_name=candidate.get("display_name", ""),
            source=candidate.get("source", ""),
            original_id=candidate.get("original_id", ""),
            predictor=candidate.get("predictor", ""),
            family=candidate.get("family", ""),
            chromosome=proposal["chromosome"],
            start=integer(proposal["start_1based"]),
            end=integer(proposal["end_1based"]),
            strand=proposal["strand"],
            coordinate_method="bowtie_exact",
            precursor=dna(proposal["b73v5_precursor_sequence"]),
            evidence_type=evidence,
            mirscore_version=candidate.get("mirscore_version", ""),
            total_reads=integer(candidate.get(f"{prefix}totReads", "")),
            supporting_libraries=integer(
                candidate.get(f"{prefix}numLibraries", candidate.get(f"{prefix}nlibraries", ""))
            ),
            precision=floating(candidate.get(f"{prefix}precision", "")),
            products=products_from_validation(candidate, evidence),
        )
        add_genomic_product_coordinates(model)
        models.append(model)
    return models


def products_by_model(path: Path) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(path):
        grouped[row["remediation_id"]].append(row)
    return grouped


def remediated_models(
    validation_path: Path,
    product_path: Path,
    branch: str,
) -> list[AcceptedModel]:
    products = products_by_model(product_path)
    models: list[AcceptedModel] = []
    for row in read_tsv(validation_path):
        status = row.get("remediation_mirscore_status", "")
        if status == "primary_mirscore_pass":
            evidence = "primary"
            prefix = "mirscore_"
        elif status == "rescue_mirscore_pass":
            evidence = "rescue"
            prefix = "mirscore_rescue_"
        else:
            continue
        remediation_id = row["remediation_id"]
        accepted_products: list[ProductEvidence] = []
        for product in products.get(remediation_id, []):
            accepted_products.append(
                ProductEvidence(
                    role=product["product_role"],
                    arm=product["arm"],
                    sequence=dna(product["mature_sequence_DNA"]),
                    relative_start=integer(product["start_1based"]),
                    relative_end=integer(product["end_1based"]),
                )
            )
        model = AcceptedModel(
            model_id=remediation_id,
            branch=branch,
            source_internal_id=row["source_internal_id"],
            display_name=row.get("display_name", ""),
            source=row.get("source", ""),
            original_id=row.get("original_id", ""),
            predictor=row.get("predictor", ""),
            family=row.get("family", ""),
            chromosome=row["chromosome"],
            start=integer(row["start_1based"]),
            end=integer(row["end_1based"]),
            strand=row["strand"],
            coordinate_method=row.get("coordinate_method", "input_remediation"),
            precursor=dna(row["b73v5_precursor_sequence"]),
            evidence_type=evidence,
            mirscore_version=row.get("mirscore_version", ""),
            total_reads=integer(row.get(f"{prefix}totReads", "")),
            supporting_libraries=integer(
                row.get(f"{prefix}numLibraries", row.get(f"{prefix}nlibraries", ""))
            ),
            precision=floating(row.get(f"{prefix}precision", "")),
            products=accepted_products,
        )
        add_genomic_product_coordinates(model)
        models.append(model)
    return models


def overlap_components(models: list[AcceptedModel]) -> list[list[AcceptedModel]]:
    axes: dict[tuple[str, str], list[AcceptedModel]] = defaultdict(list)
    for model in models:
        axes[(model.chromosome, model.strand)].append(model)
    components: list[list[AcceptedModel]] = []
    for axis in sorted(axes, key=lambda value: (natural_contig_key(value[0]), value[1])):
        ordered = sorted(
            axes[axis], key=lambda item: (item.start, item.end, item.model_id)
        )
        current: list[AcceptedModel] = []
        current_end = -1
        for model in ordered:
            if current and model.start > current_end:
                components.append(current)
                current = []
                current_end = -1
            current.append(model)
            current_end = max(current_end, model.end)
        if current:
            components.append(current)
    return sorted(
        components,
        key=lambda component: (
            natural_contig_key(component[0].chromosome),
            min(model.start for model in component),
            max(model.end for model in component),
            component[0].strand,
        ),
    )


def source_rank(source: str) -> int:
    if "expanded_mir2275_current_analysis" in source:
        return 0
    if source == "mirbase_v23":
        return 1
    if source == "mirscore_study_s4":
        return 2
    return 3


def representative_rank(model: AcceptedModel) -> tuple[object, ...]:
    coordinate_rank = {
        "bowtie_exact": 0,
        "bowtie_v3_remediation": 1,
        "minimap2_remediation": 2,
        "input_remediation": 3,
    }
    return (
        coordinate_rank.get(model.coordinate_method, 4),
        0 if model.evidence_type == "primary" else 1,
        source_rank(model.source),
        -model.total_reads,
        -model.precision,
        model.model_id,
    )


def joined(values: Iterable[str]) -> str:
    return ",".join(sorted({value for value in values if value}))


def gff_attributes(items: Iterable[tuple[str, object]]) -> str:
    return ";".join(
        f"{key}={quote(str(value), safe='_.:,|-')}"
        for key, value in items
        if value not in ("", None)
    )


def file_record_count(path: Path) -> int:
    if path.suffix in {".fa", ".fasta"}:
        with path.open("rt", encoding="utf-8") as handle:
            return sum(line.startswith(">") for line in handle)
    if path.suffix == ".gff3":
        with path.open("rt", encoding="utf-8") as handle:
            return sum(bool(line.strip()) and not line.startswith("#") for line in handle)
    if path.suffix == ".tsv":
        with path.open("rt", encoding="utf-8") as handle:
            return max(sum(1 for _ in handle) - 1, 0)
    return 1


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    candidates_list = read_tsv(args.candidate_validation)
    candidates = {row["internal_id"]: row for row in candidates_list}
    if len(candidates) != len(candidates_list):
        raise ValueError("Candidate validation contains duplicate internal IDs")

    models = exact_models(candidates, read_tsv(args.coordinate_proposals))
    models.extend(
        remediated_models(
            args.remediation_validation,
            args.remediation_products,
            "B73v5_coordinate_remediation",
        )
    )
    models.extend(
        remediated_models(
            args.input_remediation_validation,
            args.input_remediation_products,
            "input_incompatibility_remediation",
        )
    )
    if len({model.model_id for model in models}) != len(models):
        raise ValueError("Accepted model IDs are not unique")

    components = overlap_components(models)
    prefix = args.release_prefix
    locus_rows: list[dict[str, object]] = []
    evidence_rows: list[dict[str, object]] = []
    product_rows: list[dict[str, object]] = []
    gff_lines: list[str] = [
        "##gff-version 3",
        "##assembly Zm-B73-REFERENCE-NAM-5.0",
        "##annotation-source B73-v5 bona-fide MIR v1; miRScore 0.3.7",
    ]
    precursor_fasta: list[str] = []
    mature_fasta: list[str] = []
    candidate_to_loci: dict[str, set[str]] = defaultdict(set)
    global_sequences = sorted(
        {product.sequence for model in models for product in model.products}
    )
    sequence_groups = {
        sequence: f"MSEQ{index:06d}" for index, sequence in enumerate(global_sequences, 1)
    }

    for locus_number, component in enumerate(components, 1):
        locus_id = f"ZmaB73v5-MIRV1-{locus_number:06d}"
        start = min(model.start for model in component)
        end = max(model.end for model in component)
        representative = min(component, key=representative_rank)
        for transcript_number, model in enumerate(
            sorted(component, key=lambda item: (item.start, item.end, representative_rank(item))),
            1,
        ):
            model.release_locus_id = locus_id
            model.release_transcript_id = f"{locus_id}-T{transcript_number:03d}"
            candidate_to_loci[model.source_internal_id].add(locus_id)

        known_families = sorted(
            {
                model.family
                for model in component
                if model.family and model.family != "novel_or_unresolved"
            }
        )
        all_families = sorted({model.family for model in component if model.family})
        locus_family = ",".join(known_families or all_families or ["novel_or_unresolved"])
        locus_name = representative.display_name or representative.original_id or locus_id

        product_groups: dict[str, list[tuple[AcceptedModel, ProductEvidence]]] = defaultdict(list)
        for model in component:
            for product in model.products:
                product_groups[product.sequence].append((model, product))

        locus_product_records: list[dict[str, object]] = []
        for product_number, sequence in enumerate(
            sorted(
                product_groups,
                key=lambda value: (
                    min(item[1].genomic_start for item in product_groups[value]),
                    min(item[1].genomic_end for item in product_groups[value]),
                    value,
                ),
            ),
            1,
        ):
            evidence = product_groups[sequence]
            product_id = f"ZmaB73v5-miRV1-{locus_number:06d}-p{product_number:02d}"
            coordinates = sorted(
                {
                    f"{model.chromosome}:{product.genomic_start}-{product.genomic_end}:{model.strand}"
                    for model, product in evidence
                }
            )
            record = {
                "product_id": product_id,
                "release_locus_id": locus_id,
                "family": locus_family,
                "product_roles": joined(product.role for _, product in evidence),
                "arms": joined(product.arm for _, product in evidence),
                "mature_sequence_DNA": sequence,
                "mature_sequence_RNA": sequence.replace("T", "U"),
                "length_nt": len(sequence),
                "sequence_equivalence_group": sequence_groups[sequence],
                "genomic_locations": ",".join(coordinates),
                "accepted_model_ids": joined(model.model_id for model, _ in evidence),
                "source_internal_ids": joined(
                    model.source_internal_id for model, _ in evidence
                ),
                "sources": joined(model.source for model, _ in evidence),
                "original_ids": joined(model.original_id for model, _ in evidence),
                "predictors": joined(model.predictor for model, _ in evidence),
                "evidence_types": joined(model.evidence_type for model, _ in evidence),
            }
            product_rows.append(record)
            locus_product_records.append(record)
            mature_fasta.extend([f">{product_id}", sequence])

        locus_rows.append(
            {
                "release_locus_id": locus_id,
                "display_name": locus_name,
                "chromosome": representative.chromosome,
                "start_1based": start,
                "end_1based": end,
                "strand": representative.strand,
                "span_nt": end - start + 1,
                "representative_model_id": representative.model_id,
                "representative_precursor_start_1based": representative.start,
                "representative_precursor_end_1based": representative.end,
                "representative_precursor_length_nt": len(representative.precursor),
                "families": locus_family,
                "accepted_model_count": len(component),
                "accepted_product_count": len(locus_product_records),
                "source_internal_ids": joined(model.source_internal_id for model in component),
                "sources": joined(model.source for model in component),
                "original_ids": joined(model.original_id for model in component),
                "predictors": joined(model.predictor for model in component),
                "evidence_types": joined(model.evidence_type for model in component),
                "coordinate_methods": joined(model.coordinate_method for model in component),
                "representative_precursor_sequence_DNA": representative.precursor,
            }
        )
        precursor_fasta.extend([f">{locus_id}", representative.precursor])

        gene_attributes = gff_attributes(
            [
                ("ID", locus_id),
                ("Name", locus_name),
                ("gene_biotype", "miRNA"),
                ("family", locus_family),
                ("sources", joined(model.source for model in component)),
                ("release", prefix),
            ]
        )
        gff_lines.append(
            "\t".join(
                [
                    representative.chromosome,
                    prefix,
                    "gene",
                    str(start),
                    str(end),
                    ".",
                    representative.strand,
                    ".",
                    gene_attributes,
                ]
            )
        )
        sequence_to_reference_product = {
            str(row["mature_sequence_DNA"]): str(row["product_id"])
            for row in locus_product_records
        }
        for model in sorted(component, key=lambda item: item.release_transcript_id):
            transcript_attributes = gff_attributes(
                [
                    ("ID", model.release_transcript_id),
                    ("Parent", locus_id),
                    ("Name", model.display_name or model.original_id or model.model_id),
                    ("source_internal_id", model.source_internal_id),
                    ("original_id", model.original_id),
                    ("predictor", model.predictor),
                    ("miRScore_evidence", model.evidence_type),
                    ("coordinate_method", model.coordinate_method),
                    ("representative", str(model is representative).lower()),
                ]
            )
            gff_lines.append(
                "\t".join(
                    [
                        model.chromosome,
                        prefix,
                        "miRNA_primary_transcript",
                        str(model.start),
                        str(model.end),
                        ".",
                        model.strand,
                        ".",
                        transcript_attributes,
                    ]
                )
            )
            for product_number, product in enumerate(model.products, 1):
                child_id = f"{model.release_transcript_id}-P{product_number:02d}"
                child_attributes = gff_attributes(
                    [
                        ("ID", child_id),
                        ("Parent", model.release_transcript_id),
                        ("reference_product_id", sequence_to_reference_product[product.sequence]),
                        ("product_role", product.role),
                        ("arm", product.arm),
                        ("sequence_equivalence_group", sequence_groups[product.sequence]),
                    ]
                )
                gff_lines.append(
                    "\t".join(
                        [
                            model.chromosome,
                            prefix,
                            "miRNA",
                            str(product.genomic_start),
                            str(product.genomic_end),
                            ".",
                            model.strand,
                            ".",
                            child_attributes,
                        ]
                    )
                )

    for model in sorted(
        models,
        key=lambda item: (
            natural_contig_key(item.chromosome),
            item.start,
            item.end,
            item.strand,
            item.model_id,
        ),
    ):
        evidence_rows.append(
            {
                "accepted_model_id": model.model_id,
                "release_locus_id": model.release_locus_id,
                "release_transcript_id": model.release_transcript_id,
                "branch": model.branch,
                "source_internal_id": model.source_internal_id,
                "display_name": model.display_name,
                "source": model.source,
                "original_id": model.original_id,
                "predictor": model.predictor,
                "family": model.family,
                "chromosome": model.chromosome,
                "start_1based": model.start,
                "end_1based": model.end,
                "strand": model.strand,
                "coordinate_method": model.coordinate_method,
                "accepted_evidence_type": model.evidence_type,
                "mirscore_total_reads": model.total_reads,
                "mirscore_supporting_libraries": model.supporting_libraries,
                "mirscore_precision": f"{model.precision:.6g}",
                "mirscore_version": model.mirscore_version,
                "product_position_resolution": joined(
                    product.position_resolution for product in model.products
                ),
                "precursor_sequence_DNA": model.precursor,
            }
        )

    remediation_status = {
        row["source_internal_id"]: row.get("remediation_mirscore_status", "")
        for row in read_tsv(args.remediation_validation)
    }
    input_remediation_status = {
        row["source_internal_id"]: row.get("remediation_mirscore_status", "")
        for row in read_tsv(args.input_remediation_validation)
    }
    ledger_rows: list[dict[str, object]] = []
    for candidate in candidates_list:
        internal_id = candidate["internal_id"]
        loci = sorted(candidate_to_loci.get(internal_id, set()))
        original_evidence = original_evidence_type(candidate)
        if loci:
            final_status = "released_bona_fide_B73v5_MIR"
            exclusion_reason = ""
        elif internal_id in input_remediation_status:
            final_status = "excluded"
            exclusion_reason = input_remediation_status[internal_id]
        elif internal_id in remediation_status:
            final_status = "excluded"
            exclusion_reason = remediation_status[internal_id]
        elif original_evidence:
            final_status = "excluded"
            exclusion_reason = "no_accepted_passing_B73v5_model"
        elif candidate.get("mirscore_input_status") != "eligible":
            final_status = "excluded"
            exclusion_reason = "input_incompatible_not_recovered_and_passing"
        else:
            final_status = "excluded"
            exclusion_reason = "original_mirscore_fail_not_rescued"
        ledger_rows.append(
            {
                "internal_id": internal_id,
                "display_name": candidate.get("display_name", ""),
                "source": candidate.get("source", ""),
                "original_id": candidate.get("original_id", ""),
                "predictor": candidate.get("predictor", ""),
                "family": candidate.get("family", ""),
                "original_mirscore_input_status": candidate.get("mirscore_input_status", ""),
                "original_mirscore_result": candidate.get("mirscore_result", ""),
                "original_mirscore_rescue_result": candidate.get("mirscore_rescue_result", ""),
                "B73v5_remediation_status": remediation_status.get(internal_id, ""),
                "input_remediation_status": input_remediation_status.get(internal_id, ""),
                "final_release_status": final_status,
                "release_locus_ids": ",".join(loci),
                "release_locus_count": len(loci),
                "exclusion_reason": exclusion_reason,
            }
        )

    output_paths = {
        "loci": args.output_dir / f"{prefix}.loci.tsv",
        "evidence_models": args.output_dir / f"{prefix}.accepted_models.tsv",
        "products": args.output_dir / f"{prefix}.mature_products.tsv",
        "candidate_ledger": args.output_dir / f"{prefix}.candidate_ledger.tsv",
        "precursors": args.output_dir / f"{prefix}.precursors.fa",
        "mature_products": args.output_dir / f"{prefix}.mature_products.fa",
        "gff3": args.output_dir / f"{prefix}.gff3",
        "family_summary": args.output_dir / f"{prefix}.family_summary.tsv",
        "source_summary": args.output_dir / f"{prefix}.source_summary.tsv",
        "qc": args.output_dir / f"{prefix}.qc.tsv",
        "summary": args.output_dir / f"{prefix}.summary.json",
    }
    write_tsv(output_paths["loci"], locus_rows, list(locus_rows[0]))
    write_tsv(output_paths["evidence_models"], evidence_rows, list(evidence_rows[0]))
    write_tsv(output_paths["products"], product_rows, list(product_rows[0]))
    write_tsv(output_paths["candidate_ledger"], ledger_rows, list(ledger_rows[0]))
    output_paths["precursors"].write_text("\n".join(precursor_fasta) + "\n", encoding="utf-8")
    output_paths["mature_products"].write_text("\n".join(mature_fasta) + "\n", encoding="utf-8")
    output_paths["gff3"].write_text("\n".join(gff_lines) + "\n", encoding="utf-8")

    family_to_loci: dict[str, set[str]] = defaultdict(set)
    family_to_products: dict[str, set[str]] = defaultdict(set)
    for row in locus_rows:
        for family in str(row["families"]).split(","):
            family_to_loci[family].add(str(row["release_locus_id"]))
    for row in product_rows:
        for family in str(row["family"]).split(","):
            family_to_products[family].add(str(row["product_id"]))
    family_rows = [
        {
            "family": family,
            "released_loci": len(family_to_loci[family]),
            "released_products": len(family_to_products[family]),
        }
        for family in sorted(family_to_loci)
    ]
    write_tsv(output_paths["family_summary"], family_rows, list(family_rows[0]))

    source_rows = []
    for source in sorted({model.source for model in models}):
        source_models = [model for model in models if model.source == source]
        source_rows.append(
            {
                "source": source,
                "accepted_models": len(source_models),
                "released_loci_touched": len({model.release_locus_id for model in source_models}),
                "source_candidate_ids": len({model.source_internal_id for model in source_models}),
                "primary_models": sum(model.evidence_type == "primary" for model in source_models),
                "rescue_models": sum(model.evidence_type == "rescue" for model in source_models),
            }
        )
    write_tsv(output_paths["source_summary"], source_rows, list(source_rows[0]))

    non_novel_family_conflicts = []
    for row in locus_rows:
        families = {
            family
            for family in str(row["families"]).split(",")
            if family and family != "novel_or_unresolved"
        }
        if len(families) > 1:
            non_novel_family_conflicts.append(str(row["release_locus_id"]))
    expected_gff_features = len(locus_rows) + len(models) + 2 * len(models)
    qc_rows = [
        {"check": "accepted_model_ids_unique", "observed": len({m.model_id for m in models}), "expected": len(models), "status": "pass"},
        {"check": "accepted_models_assigned_to_locus", "observed": sum(bool(m.release_locus_id) for m in models), "expected": len(models), "status": "pass"},
        {"check": "candidate_ledger_complete", "observed": len(ledger_rows), "expected": len(candidates_list), "status": "pass"},
        {"check": "precursor_fasta_records", "observed": len(precursor_fasta) // 2, "expected": len(locus_rows), "status": "pass"},
        {"check": "mature_fasta_records", "observed": len(mature_fasta) // 2, "expected": len(product_rows), "status": "pass"},
        {"check": "gff_feature_records", "observed": len(gff_lines) - 3, "expected": expected_gff_features, "status": "pass"},
        {"check": "non_novel_cross_family_overlap_components", "observed": len(non_novel_family_conflicts), "expected": 0, "status": "review" if non_novel_family_conflicts else "pass"},
        {"check": "accepted_product_sequences_ACGT_19_24nt", "observed": sum(DNA_RE.fullmatch(str(r["mature_sequence_DNA"])) is not None and 19 <= int(r["length_nt"]) <= 24 for r in product_rows), "expected": len(product_rows), "status": "pass"},
    ]
    for row in qc_rows:
        if row["check"] != "non_novel_cross_family_overlap_components" and row["observed"] != row["expected"]:
            row["status"] = "fail"
    write_tsv(output_paths["qc"], qc_rows, ["check", "observed", "expected", "status"])
    if any(row["status"] == "fail" for row in qc_rows):
        raise RuntimeError("One or more release QC invariants failed")

    status_counts = Counter(str(row["final_release_status"]) for row in ledger_rows)
    summary = {
        "script": Path(__file__).name,
        "script_version": VERSION,
        "release_prefix": prefix,
        "assembly_identity": ASSEMBLY,
        "annotation_identity": ANNOTATION,
        "candidate_ledger_records": len(ledger_rows),
        "released_source_candidate_ids": sum(bool(candidate_to_loci.get(item)) for item in candidates),
        "released_loci": len(locus_rows),
        "accepted_models": len(models),
        "accepted_model_branch_counts": dict(sorted(Counter(model.branch for model in models).items())),
        "accepted_model_evidence_counts": dict(sorted(Counter(model.evidence_type for model in models).items())),
        "accepted_products_with_corrected_mirscore_offset": sum(
            product.position_resolution
            == "unique_precursor_sequence_match_corrected_miRScore_offset"
            for model in models
            for product in model.products
        ),
        "released_locus_products": len(product_rows),
        "unique_mature_product_sequences": len(global_sequences),
        "duplicate_sequence_groups_across_loci": sum(
            len({row["release_locus_id"] for row in product_rows if row["mature_sequence_DNA"] == sequence}) > 1
            for sequence in global_sequences
        ),
        "candidate_final_status_counts": dict(sorted(status_counts.items())),
        "non_novel_cross_family_overlap_loci": non_novel_family_conflicts,
        "maximum_locus_span_nt": max(int(row["span_nt"]) for row in locus_rows),
        "maximum_accepted_models_per_locus": max(int(row["accepted_model_count"]) for row in locus_rows),
        "mir2118_released_loci": len(family_to_loci.get("miR2118", set())),
        "mir2118_released_products": len(family_to_products.get("miR2118", set())),
        "mir2275_released_loci": len(family_to_loci.get("miR2275", set())),
        "mir2275_released_products": len(family_to_products.get("miR2275", set())),
        "qc_status": "review" if non_novel_family_conflicts else "pass",
    }
    output_paths["summary"].write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    command = " ".join(shlex.quote(value) for value in [Path(__file__).name, *sys.argv[1:]])
    input_paths = [
        args.candidate_validation,
        args.coordinate_proposals,
        args.remediation_validation,
        args.remediation_products,
        args.input_remediation_validation,
        args.input_remediation_products,
        args.release_policy,
        args.release_decisions,
    ]
    manifest_rows: list[dict[str, object]] = []
    for role, paths in (("input", input_paths), ("output", list(output_paths.values()))):
        for path in paths:
            manifest_rows.append(
                {
                    "role": role,
                    "path": str(path.resolve()),
                    "sha256": sha256(path),
                    "size_bytes": path.stat().st_size,
                    "record_count": file_record_count(path),
                    "assembly_identity": ASSEMBLY,
                    "annotation_identity": ANNOTATION,
                    "tool": Path(__file__).name,
                    "tool_version": VERSION,
                    "command": command,
                    "status": "current" if summary["qc_status"] == "pass" else "review_required",
                }
            )
    write_tsv(
        args.manifest_output,
        manifest_rows,
        [
            "role",
            "path",
            "sha256",
            "size_bytes",
            "record_count",
            "assembly_identity",
            "annotation_identity",
            "tool",
            "tool_version",
            "command",
            "status",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
