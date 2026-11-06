#!/usr/bin/env python3
"""Build a B73-v5 unified, publishable small-RNA-generating-locus GFF3.

This builder joins the approved structural ncRNA reference, miRBase mature
miRNA/precursor evidence, and final manually curated PHAS locus references.
It intentionally does *not* run Phasis or ShortStack.  ShortStack-exclusive
and manually merged PHAS additions are represented as locus records only,
because no Phasis phasiRNA/register products exist for them yet.  The builder
also supports the strict miRScore-qualified MIR GFF3 and explicitly separate
exploratory candidate classes used by the B73-v5 small-RNA reference v1.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import quote


DNA = set("ACGT")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_contig(contig: str) -> str:
    """Use the MaizeGDB B73-v5 `chrN` naming used by PHAS input genome."""
    return f"chr{contig}" if contig.isdigit() else contig


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", value)


def attr_value(value: object) -> str:
    return quote(str(value), safe="A-Za-z0-9.:|_-+")


def attributes(values: dict[str, object]) -> str:
    rendered: list[str] = []
    for key, value in values.items():
        if value is None or value == "":
            continue
        if key == "Parent":
            # A literal comma is the GFF3 separator for multiple parents.
            value_text = ",".join(attr_value(item) for item in str(value).split(","))
        else:
            value_text = attr_value(value)
        rendered.append(f"{key}={value_text}")
    return ";".join(rendered)


def gff_line(
    seqid: str,
    source: str,
    feature_type: str,
    start: int,
    end: int,
    strand: str,
    attrs: dict[str, object],
    score: str = ".",
    phase: str = ".",
) -> str:
    if start < 1 or end < start:
        raise ValueError(f"Invalid GFF interval {seqid}:{start}-{end}")
    return "\t".join(
        [
            seqid,
            source,
            feature_type,
            str(start),
            str(end),
            score,
            strand,
            phase,
            attributes(attrs),
        ]
    )


def parse_attributes(raw: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for item in raw.strip().split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        fields[key] = value
    return fields


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence).upper().replace("U", "T")
                header = line[1:].split(maxsplit=1)[0]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        records[header] = "".join(sequence).upper().replace("U", "T")
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def tsv_rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def parse_exact_loci(raw: str, length: int):
    for item in raw.split(";"):
        if not item:
            continue
        match = re.fullmatch(r"([^:]+):(\d+):([+-])", item)
        if not match:
            raise ValueError(f"Unrecognized exact_B73_loci item {item!r}")
        contig, start, strand = match.groups()
        start_i = int(start)
        yield normalize_contig(contig), start_i, start_i + length - 1, strand


@dataclass
class FeatureBlock:
    seqid: str
    start: int
    end: int
    lines: list[str]
    group: str


def structural_blocks(feature_decisions: Path, trnascan_gff: Path):
    decisions = list(tsv_rows(feature_decisions))
    included_trna_ids = {
        row["feature_id"]
        for row in decisions
        if row["feature_origin"] == "tRNAscan-SE" and row["decision"] == "include"
    }
    blocks: list[FeatureBlock] = []
    counts: Counter[str] = Counter()

    # tRNAscan-SE parents and their genuine exon children are transferred as a block.
    trna_records: dict[str, tuple[list[str], list[list[str]]]] = {}
    with trnascan_gff.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed tRNAscan GFF line: {raw!r}")
            attrs = parse_attributes(fields[8])
            if fields[2] == "tRNA":
                original_id = attrs.get("ID")
                if original_id in included_trna_ids:
                    trna_records[original_id] = (fields, [])
            elif fields[2] == "exon":
                parent = attrs.get("Parent")
                if parent in included_trna_ids:
                    trna_records.setdefault(parent, ([], []))[1].append(fields)
    missing = included_trna_ids - trna_records.keys()
    if missing:
        raise ValueError(f"{len(missing)} included tRNAs are absent from tRNAscan GFF")
    for original_id, (parent, children) in trna_records.items():
        if not parent:
            raise ValueError(f"Missing tRNA parent line for {original_id}")
        attrs = parse_attributes(parent[8])
        new_parent_id = f"sRNA_tRNA__{safe_id(original_id)}"
        parent_attrs = {
            "ID": new_parent_id,
            "Name": attrs.get("Name", original_id),
            "biotype": "tRNA",
            "source_feature_id": original_id,
            "isotype": attrs.get("isotype"),
            "anticodon": attrs.get("anticodon"),
            "curation_release": "B73v5_nuclear_structural_ncRNA_v1_2026-08-03",
            "evidence_tier": "canonical",
        }
        lines = [
            gff_line(
                normalize_contig(parent[0]),
                "tRNAscan-SE_2.0.13",
                "tRNA",
                int(parent[3]),
                int(parent[4]),
                parent[6],
                parent_attrs,
                parent[5],
            )
        ]
        for child_number, child in enumerate(children, start=1):
            lines.append(
                gff_line(
                    normalize_contig(child[0]),
                    "tRNAscan-SE_2.0.13",
                    "exon",
                    int(child[3]),
                    int(child[4]),
                    child[6],
                    {
                        "ID": f"{new_parent_id}.exon{child_number}",
                        "Parent": new_parent_id,
                        "biotype": "tRNA_exon",
                        "evidence_tier": "canonical",
                    },
                )
            )
            counts["tRNA_exon"] += 1
        blocks.append(
            FeatureBlock(
                normalize_contig(parent[0]), int(parent[3]), int(parent[4]), lines, "tRNA"
            )
        )
        counts["tRNA"] += 1

    # RNAcentral features are source-reference sequences placed by exact full-length
    # mapping.  rRNA entries can be reference segments rather than one full rDNA gene.
    for row in decisions:
        if row["feature_origin"] != "RNAcentral" or row["decision"] != "include":
            continue
        rna_class = row["effective_class"]
        if rna_class not in {"rRNA", "snRNA", "snoRNA"}:
            continue
        length = int(row["sequence_length"])
        for placement_number, (seqid, start, end, strand) in enumerate(
            parse_exact_loci(row["exact_B73_loci"], length), start=1
        ):
            feature_id = f"sRNA_{rna_class}__{safe_id(row['feature_id'])}__loc{placement_number}"
            attrs = {
                "ID": feature_id,
                "Name": row["feature_id"],
                "biotype": rna_class,
                "source_feature_id": row["feature_id"],
                "source_description": row["source_description"],
                "rfam_models": row["rfam_models"],
                "curation_decision": row["decision_reason"],
                "curation_release": "B73v5_nuclear_structural_ncRNA_v1_2026-08-03",
                "coordinate_evidence": "exact_full_length_B73v5_mapping",
                "feature_scope": (
                    "exact_reference_sequence_placement"
                    if rna_class == "rRNA"
                    else "exact_full_length_RNAcentral_feature_placement"
                ),
                "evidence_tier": "canonical",
            }
            blocks.append(
                FeatureBlock(
                    seqid,
                    start,
                    end,
                    [
                        gff_line(
                            seqid,
                            "RNAcentral_v27_Rfam",
                            rna_class,
                            start,
                            end,
                            strand,
                            attrs,
                        )
                    ],
                    rna_class,
                )
            )
            counts[rna_class] += 1
    return blocks, counts


def imported_mirna_blocks(path: Path):
    """Load the strict MIR release while preserving its existing hierarchy."""
    blocks: list[FeatureBlock] = []
    counts: Counter[str] = Counter()
    current: list[str] = []
    current_seqid = ""
    current_start = 0
    current_end = 0

    def flush() -> None:
        nonlocal current, current_seqid, current_start, current_end
        if current:
            blocks.append(
                FeatureBlock(
                    current_seqid,
                    current_start,
                    current_end,
                    current,
                    "miRNA",
                )
            )
        current = []
        current_seqid = ""
        current_start = 0
        current_end = 0

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed MIR GFF3 line: {raw!r}")
            feature_type = fields[2]
            if feature_type == "gene":
                flush()
                current_seqid = normalize_contig(fields[0])
                current_start = int(fields[3])
                current_end = int(fields[4])
            if not current_seqid:
                raise ValueError("Strict MIR GFF3 must begin each block with a gene feature")
            fields[0] = normalize_contig(fields[0])
            if "evidence_tier=" not in fields[8]:
                fields[8] += ";evidence_tier=canonical"
            current.append("\t".join(fields))
            counts[f"miRNA_GFF3_{feature_type}"] += 1
    flush()
    return blocks, counts


def parse_sam_mappings(sam_gz: Path, hairpins: dict[str, str]):
    mappings: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    with gzip.open(sam_gz, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("@"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"Malformed SAM record: {raw!r}")
            query, flag_text, target, position, _, cigar = fields[:6]
            flag = int(flag_text)
            if flag & 4:
                continue
            if query not in hairpins:
                raise ValueError(f"SAM query {query} absent from hairpin FASTA")
            expected_cigar = f"{len(hairpins[query])}M"
            if cigar != expected_cigar:
                raise ValueError(f"Non-full-length hairpin alignment for {query}: {cigar}")
            optional = fields[11:]
            if "NM:i:0" not in optional:
                raise ValueError(f"Non-exact hairpin alignment for {query}")
            strand = "-" if flag & 16 else "+"
            start = int(position)
            end = start + len(hairpins[query]) - 1
            mappings[query].append((normalize_contig(target), start, end, strand))
    return mappings


def mature_positions(hairpin: str, mature: str):
    """Yield 1-based positions in forward hairpin coordinates and match orientation."""
    start = 0
    while True:
        start = hairpin.find(mature, start)
        if start < 0:
            break
        yield start + 1, "forward"
        start += 1
    rc_hairpin = reverse_complement(hairpin)
    start = 0
    while True:
        start = rc_hairpin.find(mature, start)
        if start < 0:
            break
        # Convert an interval in the reverse-complemented hairpin back to the
        # forward hairpin coordinate system.
        yield len(hairpin) - (start + len(mature)) + 1, "reverse_complement"
        start += 1


def mirna_blocks(
    hairpin_fasta: Path,
    mature_fasta: Path,
    exact_sam: Path,
    ambiguous_report: Path,
    status_out: Path,
):
    hairpins = parse_fasta(hairpin_fasta)
    mature = parse_fasta(mature_fasta)
    mappings = parse_sam_mappings(exact_sam, hairpins)
    all_hairpin_ids = {
        line.split("\t", 1)[0]
        for line in ambiguous_report.read_text(encoding="utf-8").splitlines()[1:]
        if line
    }
    all_hairpin_ids |= set(hairpins)
    records: list[dict[str, object]] = []
    for precursor_id in sorted(all_hairpin_ids):
        if precursor_id not in hairpins:
            status = "not_mapped_exactly_IUPAC_ambiguous_precursor"
            records.append({"precursor_id": precursor_id, "status": status, "mapped_placements": 0})
        elif precursor_id not in mappings:
            records.append(
                {
                    "precursor_id": precursor_id,
                    "status": "not_mapped_exactly_to_B73v5",
                    "mapped_placements": 0,
                }
            )
        else:
            records.append(
                {
                    "precursor_id": precursor_id,
                    "status": "mapped_exactly_to_B73v5",
                    "mapped_placements": len(mappings[precursor_id]),
                }
            )
    status_out.parent.mkdir(parents=True, exist_ok=True)
    with status_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["precursor_id", "status", "mapped_placements"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)

    associations: dict[str, list[tuple[str, str, int, str]]] = defaultdict(list)
    for precursor_id, hairpin in hairpins.items():
        for mature_id, mature_sequence in mature.items():
            for position, orientation in mature_positions(hairpin, mature_sequence):
                associations[precursor_id].append(
                    (mature_id, mature_sequence, position, orientation)
                )
    blocks: list[FeatureBlock] = []
    counts: Counter[str] = Counter()
    for precursor_id in sorted(mappings):
        for placement_number, (seqid, start, end, strand) in enumerate(
            mappings[precursor_id], start=1
        ):
            parent_id = f"miRNA_precursor__{safe_id(precursor_id)}__loc{placement_number}"
            children: list[str] = []
            for child_number, (mature_id, sequence, hairpin_pos, orientation) in enumerate(
                associations[precursor_id], start=1
            ):
                if strand == "+":
                    child_start = start + hairpin_pos - 1
                    child_end = child_start + len(sequence) - 1
                else:
                    child_end = end - hairpin_pos + 1
                    child_start = child_end - len(sequence) + 1
                children.append(
                    gff_line(
                        seqid,
                        "miRBase_22.1",
                        "miRNA",
                        child_start,
                        child_end,
                        strand,
                        {
                            "ID": f"{parent_id}.mature{child_number}",
                            "Parent": parent_id,
                            "Name": mature_id,
                            "biotype": "miRNA",
                            "mature_sequence": sequence,
                            "mature_match_orientation_in_hairpin": orientation,
                            "hairpin_position_1based": hairpin_pos,
                        },
                    )
                )
                counts["miRNA"] += 1
            parent = gff_line(
                seqid,
                "miRBase_22.1",
                "miRNA_primary_transcript",
                start,
                end,
                strand,
                {
                    "ID": parent_id,
                    "Name": precursor_id,
                    "biotype": "miRNA_precursor",
                    "source_database": "miRBase_release_22.1",
                    "coordinate_evidence": "exact_full_length_B73v5_mapping",
                    "mature_product_count": len(children),
                },
            )
            blocks.append(FeatureBlock(seqid, start, end, [parent, *children], "miRNA_precursor"))
            counts["miRNA_primary_transcript"] += 1
    counts["miRNA_unmapped_or_ambiguous_precursor"] = sum(
        1 for row in records if row["status"] != "mapped_exactly_to_B73v5"
    )
    return blocks, counts


@dataclass
class PhasProduct:
    position: int
    sequence: str
    strand: str
    parent_windows: set[str] = field(default_factory=set)
    source_rows: int = 0
    libraries: set[str] = field(default_factory=set)
    abundance_sum: float = 0.0
    register_classes: set[str] = field(default_factory=set)
    expected_register_positions: set[int] = field(default_factory=set)


def load_phas_reference(path: Path):
    loci: dict[str, dict[str, str]] = {}
    for row in tsv_rows(path):
        if row["manual_curation_class"] != "bona_fide_PHAS":
            raise ValueError(f"Unexpected non-bona-fide row in final reference: {row}")
        loci[row["reference_locus_id"]] = row
    return loci


def phas_blocks(reference_path: Path, products_path: Path, phase_nt: int):
    loci = load_phas_reference(reference_path)
    locus_to_public = {
        locus_id: (row.get("frozen_locus_id") or locus_id) for locus_id, row in loci.items()
    }
    identifier_to_locus = {
        row["source_identifier"]: locus_id
        for locus_id, row in loci.items()
        if row["source_method"] == "Phasis"
    }
    windows: dict[tuple[str, str, str, str, str, str], dict[str, object]] = {}
    products: dict[tuple[str, str, int, str], PhasProduct] = {}
    ignored_rows = 0
    for row in tsv_rows(products_path):
        if int(row["phase"]) != phase_nt:
            raise ValueError(f"Unexpected phase in {products_path}: {row['phase']}")
        locus_id = identifier_to_locus.get(row["identifier"])
        if locus_id is None:
            ignored_rows += 1
            continue
        sequence = row["tag_seq"].upper().replace("U", "T")
        if len(sequence) != phase_nt or set(sequence) - DNA:
            raise ValueError(f"Invalid {phase_nt}-nt phasiRNA product: {sequence!r}")
        strand = "+" if row["strand"] == "w" else "-" if row["strand"] == "c" else None
        if strand is None:
            raise ValueError(f"Unexpected PHAS strand: {row['strand']!r}")
        window_key = (
            locus_id,
            strand,
            row["window_unit_id"],
            row["window_unit_role"],
            row["window_unit_rank"],
            row["window_unit_shift_nt"],
        )
        window_id = (
            f"{safe_id(locus_to_public[locus_id])}.window__{safe_id(row['window_unit_id'])}"
            f"__{safe_id(row['window_unit_role'])}__rank{safe_id(row['window_unit_rank'])}"
            f"__shift{safe_id(row['window_unit_shift_nt'])}__{strand}"
        )
        windows.setdefault(window_key, {"id": window_id, "products": set()})
        product_key = (locus_id, strand, int(row["observed_pos"]), sequence)
        product = products.setdefault(
            product_key,
            PhasProduct(int(row["observed_pos"]), sequence, strand),
        )
        product.parent_windows.add(window_id)
        product.source_rows += 1
        product.libraries.add(row["alib"])
        product.abundance_sum += float(row["abun"])
        product.register_classes.add(row["register_class"])
        product.expected_register_positions.add(int(row["expected_register_pos"]))
        windows[window_key]["products"].add(product_key)

    products_by_locus: dict[str, list[PhasProduct]] = defaultdict(list)
    product_keys_by_locus: dict[str, list[tuple[str, str, int, str]]] = defaultdict(list)
    for product_key, product in products.items():
        products_by_locus[product_key[0]].append(product)
        product_keys_by_locus[product_key[0]].append(product_key)
    windows_by_locus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for window_key, window in windows.items():
        windows_by_locus[window_key[0]].append(window)

    blocks: list[FeatureBlock] = []
    counts: Counter[str] = Counter()
    source_method_counts: Counter[str] = Counter()
    expanded_parent_count = 0
    for locus_id, row in loci.items():
        public_locus_id = row.get("frozen_locus_id") or locus_id
        seqid = row["chromosome"]
        curated_start = int(row["start_1based"])
        curated_end = int(row["end_1based"])
        locus_products = products_by_locus[locus_id]
        span_start = min([curated_start, *[p.position for p in locus_products]])
        span_end = max([curated_end, *[p.position + phase_nt - 1 for p in locus_products]])
        adjusted = span_start != curated_start or span_end != curated_end
        if adjusted:
            expanded_parent_count += 1
        phas_id = f"sRNA_{phase_nt}PHAS__{safe_id(public_locus_id)}"
        parent_attrs = {
            "ID": phas_id,
            "Name": public_locus_id,
            "biotype": f"{phase_nt}-PHAS",
            "evidence_tier": "canonical",
            "phase_length_nt": phase_nt,
            "source_reference_locus_id": locus_id,
            "source_method": row["source_method"],
            "source_identifier": row["source_identifier"],
            "manual_curation_class": row["manual_curation_class"],
            "manual_curation_basis": row["manual_curation_basis"],
            "phasis_evidence_category": row["phasis_evidence_category"],
            "shortstack_phase_score": row["shortstack_phase_score"],
            "shortstack_reads": row["shortstack_reads"],
            "maximum_Peak_Howell_score": row.get("maximum_Peak_Howell_score"),
            "Howell_support_call_count": row.get("Howell_support_call_count"),
            "freeze_release": row.get("freeze_release"),
            "freeze_status": row.get("freeze_status"),
            "curated_start_1based": curated_start,
            "curated_end_1based": curated_end,
            "gff_span_expanded_for_phasiRNA_children": str(adjusted).lower(),
            "phasiRNA_product_status": (
                "existing_Phasis_combined_product_table" if locus_products else "not_yet_defined"
            ),
        }
        lines = [
            gff_line(
                seqid,
                "manual_Phasis_ShortStack_curation",
                "ncRNA_gene",
                span_start,
                span_end,
                ".",
                parent_attrs,
            )
        ]
        for window in sorted(windows_by_locus[locus_id], key=lambda item: str(item["id"])):
            keys = list(window["products"])
            window_products = [products[key] for key in keys]
            window_start = min(p.position for p in window_products)
            window_end = max(p.position + phase_nt - 1 for p in window_products)
            window_id = str(window["id"])
            lines.append(
                gff_line(
                    seqid,
                    "Phasis_2.8",
                    "transcript",
                    window_start,
                    window_end,
                    ".",
                    {
                        "ID": window_id,
                        "Parent": phas_id,
                        "biotype": "phasing_window",
                        "evidence_tier": "canonical",
                        "phase_length_nt": phase_nt,
                        "window_unit_id": window_id.split(".window__", 1)[1].split("__", 1)[0],
                        "window_definition": "alternative_or_nonoverlapping_phasing_window",
                    },
                )
            )
            counts[f"{phase_nt}-phasing_window"] += 1
        for product_key in sorted(
            product_keys_by_locus[locus_id], key=lambda key: (key[2], key[1], key[3])
        ):
            product = products[product_key]
            product_id = (
                f"{safe_id(public_locus_id)}.phasiRNA__{product.strand}__{product.position}"
                f"__{product.sequence}"
            )
            lines.append(
                gff_line(
                    seqid,
                    "Phasis_2.8",
                    "exon",
                    product.position,
                    product.position + phase_nt - 1,
                    product.strand,
                    {
                        "ID": product_id,
                        "Parent": ",".join(sorted(product.parent_windows)),
                        "Name": product.sequence,
                        "biotype": "phasiRNA",
                        "evidence_tier": "canonical",
                        "phase_length_nt": phase_nt,
                        "register_class": ",".join(sorted(product.register_classes)),
                        "expected_register_positions": ",".join(
                            str(value) for value in sorted(product.expected_register_positions)
                        ),
                        "source_product_rows": product.source_rows,
                        "source_library_count": len(product.libraries),
                        "source_abundance_sum": f"{product.abundance_sum:.12g}",
                    },
                )
            )
            counts[f"{phase_nt}-phasiRNA"] += 1
        blocks.append(FeatureBlock(seqid, span_start, span_end, lines, f"{phase_nt}-PHAS"))
        counts[f"{phase_nt}-PHAS"] += 1
        source_method_counts[row["source_method"]] += 1
    counts["ignored_unapproved_product_rows"] = ignored_rows
    counts["parents_with_expanded_gff_span"] = expanded_parent_count
    for method, count in source_method_counts.items():
        counts[f"{phase_nt}-PHAS_{method}"] = count
    return blocks, counts


def unsupported_mirna_blocks(path: Path):
    blocks: list[FeatureBlock] = []
    counts: Counter[str] = Counter()
    for row in tsv_rows(path):
        model_id = row["coordinate_model_id"]
        parent_id = f"sRNA_unsupported_MIR__{safe_id(model_id)}"
        seqid = normalize_contig(row["chromosome"])
        start = int(row["start_1based"])
        end = int(row["end_1based"])
        strand = row["strand"]
        gene = gff_line(
            seqid,
            "B73v5_MIR_candidate_review",
            "ncRNA_gene",
            start,
            end,
            strand,
            {
                "ID": parent_id,
                "Name": row["display_name"],
                "biotype": "unsupported_MIR_candidate",
                "evidence_tier": "exploratory",
                "family": row["family"],
                "source_internal_id": row["internal_id"],
                "coordinate_method": row["coordinate_method"],
                "coordinate_status": row["coordinate_status"],
                "exclusion_reason": row["exclusion_reason"],
            },
        )
        transcript = gff_line(
            seqid,
            "B73v5_MIR_candidate_review",
            "ncRNA",
            start,
            end,
            strand,
            {
                "ID": f"{parent_id}.candidate_precursor",
                "Parent": parent_id,
                "biotype": "unsupported_MIR_candidate_precursor_model",
                "evidence_tier": "exploratory",
            },
        )
        blocks.append(FeatureBlock(seqid, start, end, [gene, transcript], "unsupported_MIR"))
        counts["unsupported_MIR_candidate_coordinate_model"] += 1
    return blocks, counts


def phas_like_blocks(path: Path, phase_nt: int):
    blocks: list[FeatureBlock] = []
    counts: Counter[str] = Counter()
    for row in tsv_rows(path):
        locus_id = row["phas_like_locus_id"]
        seqid = normalize_contig(row["chromosome"])
        start = int(row["start_1based"])
        end = int(row["end_1based"])
        attrs = {
            "ID": f"sRNA_{phase_nt}PHASLIKE__{safe_id(locus_id)}",
            "Name": locus_id,
            "biotype": f"{phase_nt}-PHAS-like",
            "evidence_tier": "exploratory",
            "phase_length_nt": phase_nt,
            "source_identifier": row.get("source_identifier"),
            "evidence_modes": row.get("evidence_modes"),
            "maximum_Peak_Howell_score": row.get("maximum_Peak_Howell_score"),
            "manual_curation_basis": row.get("manual_curation_basis"),
            "manual_curation_status": row.get("manual_curation_status"),
            "reference_status": row.get("reference_status")
            or row.get("exploratory_class_status"),
            "product_semantics": row.get("product_semantics"),
        }
        line = gff_line(
            seqid,
            "Phasis_candidate_review",
            "ncRNA_gene",
            start,
            end,
            ".",
            attrs,
        )
        blocks.append(FeatureBlock(seqid, start, end, [line], f"{phase_nt}-PHAS-like"))
        counts[f"{phase_nt}-PHAS-like"] += 1
    return blocks, counts


def write_summary(
    path: Path,
    counts: Counter[str],
    inputs: list[Path],
    gff_path: Path,
    feature_lines: int,
):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"gff3_path\t{gff_path}\n")
        handle.write(f"gff3_sha256\t{sha256(gff_path)}\n")
        handle.write(f"gff3_feature_lines\t{feature_lines}\n")
        for key in sorted(counts):
            handle.write(f"{key}\t{counts[key]}\n")
        for input_path in inputs:
            handle.write(f"input_sha256:{input_path.name}\t{sha256(input_path)}\n")


def write_manifest(path: Path, output: Path, inputs: list[Path], release_label: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "record_type\tpath\tsha256\tassembly_identity\tannotation_identity\t"
            "command_or_source\tstatus\n"
        )
        handle.write(
            f"unified_gff3\t{output}\t{sha256(output)}\t"
            "Zm-B73-REFERENCE-NAM-5.0; B73 RefGen v5\t"
            "strict miRScore-qualified MIR v1; approved structural-ncRNA v1; "
            "frozen PHAS v1; exploratory MIR and PHAS-like candidates\t"
            f"scripts/build_b73v5_unified_srna_gff3.py ({release_label})\tcomplete\n"
        )
        for input_path in inputs:
            handle.write(
                f"input\t{input_path}\t{sha256(input_path)}\t"
                "Zm-B73-REFERENCE-NAM-5.0; B73 RefGen v5 or source sequence\t"
                "see release README\tbuild input\tverified_at_build\n"
            )


def write_gff3(path: Path, blocks: list[FeatureBlock], release_label: str) -> int:
    ordered = sorted(blocks, key=lambda block: (block.seqid, block.start, block.end, block.group))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write("##species https://identifiers.org/taxonomy:4577\n")
        handle.write("##assembly Zm-B73-REFERENCE-NAM-5.0\n")
        handle.write(f"##annotation-version {release_label}\n")
        handle.write(
            "# PHAS transcript children are Phasis phasing windows; phasiRNA exon "
            "children can have multiple Parent windows.\n"
        )
        handle.write(
            "# ShortStack/manual-only PHAS loci have no inferred phasing-window or "
            "phasiRNA children.\n"
        )
        for block in ordered:
            for line in block.lines:
                handle.write(line)
                handle.write("\n")
    return sum(len(block.lines) for block in ordered)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structural-feature-decisions", type=Path, required=True)
    parser.add_argument("--trnascan-gff3", type=Path, required=True)
    parser.add_argument("--mirna-gff3", type=Path)
    parser.add_argument("--mirna-hairpin-fasta", type=Path)
    parser.add_argument("--mirna-mature-fasta", type=Path)
    parser.add_argument("--mirna-exact-sam", type=Path)
    parser.add_argument("--mirna-ambiguous-report", type=Path)
    parser.add_argument("--unsupported-mirna-coordinate-models", type=Path)
    parser.add_argument("--21-phas-like-loci", dest="phas_like_21", type=Path)
    parser.add_argument("--24-phas-like-loci", dest="phas_like_24", type=Path)
    parser.add_argument("--21-reference", dest="reference_21", type=Path, required=True)
    parser.add_argument("--24-reference", dest="reference_24", type=Path, required=True)
    parser.add_argument("--21-products", dest="products_21", type=Path, required=True)
    parser.add_argument("--24-products", dest="products_24", type=Path, required=True)
    parser.add_argument("--output-gff3", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mirna-coordinate-status", type=Path, required=True)
    parser.add_argument("--output-mirna-gff3", type=Path)
    parser.add_argument("--output-structural-gff3", type=Path)
    parser.add_argument("--output-21-phas-gff3", type=Path)
    parser.add_argument("--output-24-phas-gff3", type=Path)
    parser.add_argument("--output-exploratory-gff3", type=Path)
    parser.add_argument("--release-label", required=True)
    args = parser.parse_args()

    structural, structural_counts = structural_blocks(
        args.structural_feature_decisions, args.trnascan_gff3
    )
    if args.mirna_gff3:
        mirna, mirna_counts = imported_mirna_blocks(args.mirna_gff3)
        args.mirna_coordinate_status.parent.mkdir(parents=True, exist_ok=True)
        args.mirna_coordinate_status.write_text(
            "status\tdetail\nnot_applicable\tstrict_MIR_GFF3_imported_directly\n",
            encoding="utf-8",
        )
    else:
        legacy_mirna_inputs = [
            args.mirna_hairpin_fasta,
            args.mirna_mature_fasta,
            args.mirna_exact_sam,
            args.mirna_ambiguous_report,
        ]
        if any(value is None for value in legacy_mirna_inputs):
            parser.error(
                "provide --mirna-gff3 or all legacy MIR FASTA/SAM/report arguments"
            )
        mirna, mirna_counts = mirna_blocks(
            args.mirna_hairpin_fasta,
            args.mirna_mature_fasta,
            args.mirna_exact_sam,
            args.mirna_ambiguous_report,
            args.mirna_coordinate_status,
        )
    phas21, phas21_counts = phas_blocks(args.reference_21, args.products_21, 21)
    phas24, phas24_counts = phas_blocks(args.reference_24, args.products_24, 24)
    exploratory: list[FeatureBlock] = []
    exploratory_counts: Counter[str] = Counter()
    if args.unsupported_mirna_coordinate_models:
        blocks, counts_part = unsupported_mirna_blocks(args.unsupported_mirna_coordinate_models)
        exploratory.extend(blocks)
        exploratory_counts.update(counts_part)
    if args.phas_like_21:
        blocks, counts_part = phas_like_blocks(args.phas_like_21, 21)
        exploratory.extend(blocks)
        exploratory_counts.update(counts_part)
    if args.phas_like_24:
        blocks, counts_part = phas_like_blocks(args.phas_like_24, 24)
        exploratory.extend(blocks)
        exploratory_counts.update(counts_part)
    all_blocks = [*structural, *mirna, *phas21, *phas24, *exploratory]
    feature_lines = write_gff3(args.output_gff3, all_blocks, args.release_label)
    for path, blocks in [
        (args.output_mirna_gff3, mirna),
        (args.output_structural_gff3, structural),
        (args.output_21_phas_gff3, phas21),
        (args.output_24_phas_gff3, phas24),
        (args.output_exploratory_gff3, exploratory),
    ]:
        if path:
            write_gff3(path, blocks, args.release_label)
    counts = Counter()
    counts.update(structural_counts)
    counts.update(mirna_counts)
    counts.update(phas21_counts)
    counts.update(phas24_counts)
    counts.update(exploratory_counts)
    inputs = [
        args.structural_feature_decisions,
        args.trnascan_gff3,
        args.reference_21,
        args.reference_24,
        args.products_21,
        args.products_24,
    ]
    if args.mirna_gff3:
        inputs.append(args.mirna_gff3)
    else:
        inputs.extend(
            [
                args.mirna_hairpin_fasta,
                args.mirna_mature_fasta,
                args.mirna_exact_sam,
                args.mirna_ambiguous_report,
            ]
        )
    inputs.extend(
        value
        for value in [
            args.unsupported_mirna_coordinate_models,
            args.phas_like_21,
            args.phas_like_24,
        ]
        if value is not None
    )
    write_summary(args.summary, counts, inputs, args.output_gff3, feature_lines)
    write_manifest(args.manifest, args.output_gff3, inputs, args.release_label)
    print(f"output_gff3={args.output_gff3}")
    print(f"feature_blocks={len(all_blocks)}")
    print(f"feature_lines={feature_lines}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
