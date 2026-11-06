#!/usr/bin/env python3
"""Feature-ID-aware validation of the B73v5 curated sRNA reference (Version 1)
against a local, faidx-indexed B73v5 genome FASTA.

NEW VALIDATION CODE written 2026-11-05. This is NOT the original annotation-
building code (that code is not in the release; see provenance/
annotation_build_manifest.tsv -> scripts/build_b73v5_unified_srna_gff3.py).

Read-only with respect to every input. Requires only the Python 3 standard
library and a .fai index next to the genome FASTA (samtools faidx).

For every coordinate-bearing feature in the combined GFF3 it:
  1. checks the interval lies inside the genome contig (coordinate compatibility);
  2. resolves the feature's OWN expected sequence through release identifiers
     (GFF Name / reference_product_id / precursor ID / structural membership
     table), never through class-level sequence membership;
  3. compares that sequence with the genome at the feature's coordinates in the
     expected orientation;
  4. checks auxiliary links (FASTA record, sequence_equivalence membership,
     Phasis observation rows, listed genomic placements, parent containment).

Sequence statuses:
  exact_for_ID                 feature's own sequence == genome, expected orientation
  reverse_complement_for_ID    feature's own sequence == genome, opposite orientation only
  only_other_ID_at_coordinate  genome sequence matches a different record of the class FASTA
  absent                       no match at all
  not_evaluable                no feature-specific sequence exists in the release

Orientation conventions (documented, then tested by the data):
  PHAS phasiRNA exons  : expected sequence = GFF Name, REFERENCE-FORWARD for both
                         strands (Phasis stores the SAM SEQ field verbatim; see report).
  miRNA, structural ncRNA: expected sequence = transcript orientation (GFF strand).
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import os
import platform
import subprocess
import sys
from urllib.parse import unquote

csv.field_size_limit(sys.maxsize)
RC_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(seq):
    return seq.translate(RC_TABLE)[::-1]


def sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')


def read_tsv(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def read_fasta(path):
    records, order, name, parts = {}, [], None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(parts).upper()
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f'duplicate FASTA id {name} in {path}')
                order.append(name)
                parts = []
            elif line:
                parts.append(line.strip())
        if name is not None:
            records[name] = ''.join(parts).upper()
    return records


class IndexedFasta:
    """Minimal random-access reader using a samtools .fai index."""

    def __init__(self, path):
        self.path = path
        self.index = {}
        with open(path + '.fai') as handle:
            for line in handle:
                name, length, offset, linebases, linewidth = line.split('\t')[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self.handle = open(path, 'rb')

    def length(self, chrom):
        return self.index[chrom][0] if chrom in self.index else None

    def fetch(self, chrom, start1, end1):
        length, offset, lb, lw = self.index[chrom]
        if start1 < 1 or end1 > length or end1 < start1:
            raise ValueError(f'{chrom}:{start1}-{end1} outside contig length {length}')
        s0, e0 = start1 - 1, end1 - 1
        first = offset + (s0 // lb) * lw + s0 % lb
        last = offset + (e0 // lb) * lw + e0 % lb
        self.handle.seek(first)
        raw = self.handle.read(last - first + 1).decode('ascii')
        seq = raw.replace('\n', '').replace('\r', '').upper()
        if len(seq) != end1 - start1 + 1:
            raise ValueError(f'extraction length mismatch at {chrom}:{start1}-{end1}')
        return seq


def parse_attrs(text):
    out = {}
    for field in text.strip().split(';'):
        if not field:
            continue
        key, _, value = field.partition('=')
        out[key] = unquote(value)
    return out


def read_gff(path):
    features = []
    with open_text(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9:
                raise ValueError(f'{path}:{lineno}: expected 9 columns')
            features.append({
                'lineno': lineno, 'seqid': f[0], 'source': f[1], 'type': f[2],
                'start': int(f[3]), 'end': int(f[4]), 'strand': f[6],
                'attrs': parse_attrs(f[8]),
            })
    return features


def read_sam_placements(path):
    """query -> list of (rname, pos1, reverse_flag). Read-only use of existing SAM."""
    placements = collections.defaultdict(list)
    sq = set()
    if not path:
        return None, sq
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):
                if line.startswith('@SQ'):
                    for tok in line.split('\t'):
                        if tok.startswith('SN:'):
                            sq.add(tok[3:].strip())
                continue
            f = line.split('\t', 4)
            flag = int(f[1])
            if flag & 4:
                continue
            placements[f[0]].append((f[2], int(f[3]), bool(flag & 16)))
    return placements, sq


def norm_chrom(c, genome):
    if genome.length(c) is not None:
        return c
    if genome.length('chr' + c) is not None:
        return 'chr' + c
    return c


def parse_locus_list(text, genome):
    """'6:16754812:-;1:5-90:+' or 'chr2:6549441-6549461:-,..' -> set of (chrom, start, end|None, strand)."""
    out = set()
    if not text:
        return out
    for item in text.replace(',', ';').split(';'):
        item = item.strip()
        if not item:
            continue
        head, _, strand = item.rpartition(':')
        chrom, _, coords = head.rpartition(':')
        chrom = norm_chrom(chrom, genome)
        if '-' in coords:
            s, e = coords.split('-')
            out.add((chrom, int(s), int(e), strand))
        else:
            out.add((chrom, int(coords), None, strand))
    return out


def orientation_match(expected, fwd):
    """Compare a feature's expected sequence with the genome reference-forward sequence."""
    f_ok, r_ok = expected == fwd, expected == rc(fwd)
    if f_ok and r_ok:
        return 'both_orientations'
    if f_ok:
        return 'reference_forward'
    if r_ok:
        return 'reverse_complement'
    return 'none'


def sequence_status(match, expected_orientation, strand, fwd, other_class_seqs):
    """Map an orientation match to the handoff status vocabulary."""
    if expected_orientation == 'reference_forward':
        want = 'reference_forward'
    else:  # transcript orientation
        want = 'reverse_complement' if strand == '-' else 'reference_forward'
    if match == 'both_orientations' or match == want:
        return 'exact_for_ID'
    if match != 'none':
        return 'reverse_complement_for_ID'
    if other_class_seqs is not None and (fwd in other_class_seqs or rc(fwd) in other_class_seqs):
        return 'only_other_ID_at_coordinate'
    return 'absent'


ROW_FIELDS = [
    'class', 'evidence_tier', 'feature_type', 'feature_id', 'seqid', 'start_1based', 'end_1based',
    'strand', 'length', 'coordinate_status', 'expected_sequence_source', 'expected_sequence_id',
    'expected_orientation', 'genome_match_orientation', 'sequence_status', 'fasta_id',
    'fasta_link', 'equivalence_link', 'observation_link', 'observation_rows',
    'observation_counts_agree', 'listed_placements', 'placement_includes_feature',
    'multimapper', 'parent_containment', 'feature_pass', 'notes',
]


def new_row(feat, cls, tier):
    return {k: '' for k in ROW_FIELDS} | {
        'class': cls, 'evidence_tier': tier, 'feature_type': feat['type'],
        'feature_id': feat['attrs'].get('ID', f'line{feat["lineno"]}'),
        'seqid': feat['seqid'], 'start_1based': feat['start'], 'end_1based': feat['end'],
        'strand': feat['strand'], 'length': feat['end'] - feat['start'] + 1,
    }


def coordinate_status(feat, genome):
    length = genome.length(feat['seqid'])
    if length is None:
        return 'contig_absent_from_genome'
    if feat['start'] < 1 or feat['end'] < feat['start'] or feat['end'] > length:
        return 'outside_contig'
    return 'valid'


def bool_str(value):
    if value is None or value == '':
        return 'NA'
    return 'yes' if value else 'no'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--release', required=True, help='B73v5_curated_sRNA_reference_v1_2026-08-27 directory')
    ap.add_argument('--genome', required=True, help='uncompressed genome FASTA with .fai')
    ap.add_argument('--sam21', help='existing Bowtie1 SAM of 21_PHAS.phasiRNAs.fa vs this genome (read-only)')
    ap.add_argument('--sam24', help='existing Bowtie1 SAM of 24_PHAS.phasiRNAs.fa vs this genome (read-only)')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--skip-genome-sha', action='store_true', help='skip hashing the multi-GB genome')
    args = ap.parse_args()

    R = args.release
    os.makedirs(args.outdir, exist_ok=True)
    gpath = lambda *p: os.path.join(R, *p)
    genome = IndexedFasta(args.genome)

    gff_path = gpath('annotation', 'B73v5_sRNA_reference_v1.combined.gff3')
    feats = read_gff(gff_path)
    by_id = {f['attrs']['ID']: f for f in feats if 'ID' in f['attrs']}

    fasta = {name: read_fasta(gpath('fasta', name)) for name in [
        '21_PHAS.phasiRNAs.fa', '24_PHAS.phasiRNAs.fa', 'miRNA_bona_fide.mature_and_star.fa',
        'miRNA_bona_fide.precursors.fa', 'tRNA.genomic.fa', 'tRNA.CCA_aware_quantification.fa',
        'rRNA.fa', 'snRNA.fa', 'snoRNA.fa', 'miRNA_unsupported_candidates.fa',
        '21_PHAS_like.candidate_products.fa', '24_PHAS_like.candidate_products.fa']}
    seq_to_ids = {k: collections.defaultdict(list) for k in fasta}
    for k, recs in fasta.items():
        for rid, s in recs.items():
            seq_to_ids[k][s].append(rid)
    # full class sequence sets; when a feature's own sequence matches in neither orientation,
    # membership in the full set necessarily means another record of the class FASTA
    seqsets = {k: set(v.values()) for k, v in fasta.items()}

    seqeq = read_tsv(gpath('tables', 'sequence_equivalence.tsv'))
    seqeq_by_id = {r['sequence_reference_id']: r for r in seqeq}
    membership = read_tsv(gpath('tables', 'structural_ncRNA.sequence_membership.tsv'))
    member_by_id = {r['member_id']: r for r in membership}
    mir_products = {r['product_id']: r for r in read_tsv(gpath('tables', 'miRNA_bona_fide.products.tsv'))}
    mir_loci = {r['release_locus_id']: r for r in read_tsv(gpath('tables', 'miRNA_bona_fide.loci.tsv'))}

    # Phasis observation rows: key (identifier, w/c, observed_pos, tag_seq)
    obs = {}
    obs_any = {}
    for phase in ('21', '24'):
        idx = collections.defaultdict(lambda: {'rows': 0, 'libs': set(), 'abun': 0.0})
        idx_any = collections.defaultdict(set)
        for r in read_tsv(gpath('tables', f'{phase}_PHAS.phasiRNA_observations.tsv.gz')):
            key = (r['identifier'], r['strand'], int(r['observed_pos']), r['tag_seq'].upper())
            entry = idx[key]
            entry['rows'] += 1
            entry['libs'].add(r['alib'])
            entry['abun'] += float(r['abun'])
            idx_any[key[1:]].add(r['identifier'])
        obs[phase], obs_any[phase] = idx, idx_any

    sam = {}
    for phase, path in (('21', args.sam21), ('24', args.sam24)):
        placements, sq = read_sam_placements(path)
        sam[phase] = placements
        if placements is not None:
            missing = sorted(set(genome.index) ^ sq)
            if missing:
                print(f'WARNING: SAM {phase} @SQ names differ from genome .fai: {missing[:5]}', file=sys.stderr)

    rows = []
    notes_global = []

    # ------------------------------------------------------------------ PHAS
    locus_by_name = {}
    for f in feats:
        if f['type'] == 'ncRNA_gene' and f['attrs'].get('biotype') in ('21-PHAS', '24-PHAS'):
            locus_by_name[f['attrs']['Name']] = f
    children = collections.Counter()
    for f in feats:
        for p in f['attrs'].get('Parent', '').split(','):
            if p:
                children[p] += 1

    for f in feats:
        a = f['attrs']
        if f['type'] != 'exon' or a.get('biotype') != 'phasiRNA':
            continue
        phase = a['phase_length_nt']
        cls = f'{phase}-PHAS_phasiRNA'
        row = new_row(f, cls, a.get('evidence_tier', ''))
        row['coordinate_status'] = coordinate_status(f, genome)
        notes = []
        fid = a['ID']
        locus_name, _, tail = fid.partition('.phasiRNA__')
        id_strand, id_pos, id_seq = tail.split('__')
        expected = a.get('Name', '').upper()
        if id_seq.upper() != expected:
            notes.append('ID_sequence_differs_from_Name')
        if id_strand != f['strand'] or int(id_pos) != f['start']:
            notes.append('ID_strand_or_position_differs_from_GFF_columns')
        if len(expected) != row['length']:
            notes.append('Name_length_differs_from_interval')
        row['expected_sequence_source'] = 'GFF Name (= sequence embedded in ID)'
        row['expected_sequence_id'] = expected
        row['expected_orientation'] = 'reference_forward'
        fasta_key = f'{phase}_PHAS.phasiRNAs.fa'
        if row['coordinate_status'] == 'valid':
            fwd = genome.fetch(f['seqid'], f['start'], f['end'])
            row['genome_match_orientation'] = orientation_match(expected, fwd)
            row['sequence_status'] = sequence_status(row['genome_match_orientation'], 'reference_forward',
                                                     f['strand'], fwd, seqsets[fasta_key])
        else:
            row['sequence_status'] = 'not_evaluable'
        ids = seq_to_ids[fasta_key].get(expected, [])
        row['fasta_id'] = ','.join(ids)
        row['fasta_link'] = bool_str(len(ids) == 1)
        eq = seqeq_by_id.get(ids[0]) if ids else None
        row['equivalence_link'] = bool_str(eq is not None and eq['sequence_DNA'].upper() == expected
                                           and locus_name in eq['member_ids'].split(','))
        locus = locus_by_name.get(locus_name)
        if locus is None:
            notes.append('parent_locus_not_found')
            row['parent_containment'] = 'no'
        else:
            row['parent_containment'] = bool_str(locus['seqid'] == f['seqid'] and locus['start'] <= f['start']
                                                 and f['end'] <= locus['end'])
            ident = locus['attrs'].get('source_identifier', '')
            wc = 'w' if f['strand'] == '+' else 'c'
            entry = obs[phase].get((ident, wc, f['start'], expected))
            if entry:
                row['observation_link'] = 'yes'
                row['observation_rows'] = entry['rows']
                agree = (str(entry['rows']) == a.get('source_product_rows')
                         and str(len(entry['libs'])) == a.get('source_library_count')
                         and abs(entry['abun'] - float(a.get('source_abundance_sum', 'nan'))) < 1e-6)
                row['observation_counts_agree'] = bool_str(agree)
            else:
                other = obs_any[phase].get((wc, f['start'], expected))
                row['observation_link'] = 'other_identifier' if other else 'no'
                if other:
                    notes.append('observation_under_identifier=' + '|'.join(sorted(other)))
        placements = sam[phase].get(ids[0], []) if (sam[phase] is not None and ids) else None
        if placements is not None:
            row['listed_placements'] = len(placements)
            row['placement_includes_feature'] = bool_str((f['seqid'], f['start'], False) in placements)
            row['multimapper'] = bool_str(len(placements) > 1)
        row['feature_pass'] = bool_str(
            row['coordinate_status'] == 'valid' and row['sequence_status'] == 'exact_for_ID'
            and row['fasta_link'] == 'yes' and row['equivalence_link'] == 'yes'
            and row['observation_link'] == 'yes' and row['parent_containment'] == 'yes'
            and row['placement_includes_feature'] in ('yes', ''))
        row['notes'] = ';'.join(notes)
        rows.append(row)

    # PHAS loci and windows: coordinate / structure only
    for f in feats:
        a = f['attrs']
        bt = a.get('biotype', '')
        if f['type'] == 'ncRNA_gene' and bt in ('21-PHAS', '24-PHAS', '21-PHAS-like', '24-PHAS-like'):
            locus_only = a.get('phasiRNA_product_status') == 'not_yet_defined'
            cls = bt + ('_locus_only' if locus_only else '_locus')
            row = new_row(f, cls, a.get('evidence_tier', ''))
            row['coordinate_status'] = coordinate_status(f, genome)
            row['sequence_status'] = 'not_evaluable'
            row['expected_sequence_source'] = 'none (locus interval; no locus sequence in release)'
            n_child = children.get(a['ID'], 0)
            notes = [f'direct_children={n_child}', f'source_method={a.get("source_method", "")}']
            cs, ce = a.get('curated_start_1based'), a.get('curated_end_1based')
            expanded = a.get('gff_span_expanded_for_phasiRNA_children') == 'true'
            if cs and ce and (int(cs), int(ce)) != (f['start'], f['end']):
                notes.append('GFF_span_differs_from_curated_span' + ('(declared_expanded)' if expanded else '(UNDECLARED)'))
            structure_ok = (n_child == 0) if (locus_only or bt.endswith('like')) else (n_child > 0)
            row['parent_containment'] = 'NA'
            row['feature_pass'] = bool_str(row['coordinate_status'] == 'valid' and structure_ok and
                                           'UNDECLARED' not in ''.join(notes))
            if not structure_ok:
                notes.append('unexpected_child_count')
            row['notes'] = ';'.join(notes)
            rows.append(row)
        elif f['type'] == 'transcript' and bt == 'phasing_window':
            row = new_row(f, f'{a["phase_length_nt"]}-PHAS_phasing_window', a.get('evidence_tier', ''))
            row['coordinate_status'] = coordinate_status(f, genome)
            row['sequence_status'] = 'not_evaluable'
            parent = by_id.get(a.get('Parent', ''))
            row['parent_containment'] = bool_str(parent is not None and parent['start'] <= f['start'] and f['end'] <= parent['end'])
            row['feature_pass'] = bool_str(row['coordinate_status'] == 'valid' and row['parent_containment'] == 'yes')
            rows.append(row)

    # ------------------------------------------------------------------ miRNA
    prec = fasta['miRNA_bona_fide.precursors.fa']
    mature = fasta['miRNA_bona_fide.mature_and_star.fa']
    for f in feats:
        a = f['attrs']
        if f['source'] != 'Zm_B73v5_bona_fide_miRNA_v1':
            continue
        row = new_row(f, 'miRNA_' + f['type'], a.get('evidence_tier', ''))
        row['coordinate_status'] = coordinate_status(f, genome)
        row['expected_orientation'] = 'transcript (GFF strand)'
        notes = []
        expected, region = None, (f['seqid'], f['start'], f['end'], f['strand'])
        if f['type'] == 'gene':
            locus = mir_loci.get(a['ID'])
            expected = prec.get(a['ID'])
            row['expected_sequence_source'] = 'miRNA_bona_fide.precursors.fa (same ID); compared at representative precursor coordinates'
            row['expected_sequence_id'] = a['ID']
            row['fasta_id'] = a['ID'] if expected else ''
            row['fasta_link'] = bool_str(expected is not None)
            if locus:
                region = (locus['chromosome'], int(locus['representative_precursor_start_1based']),
                          int(locus['representative_precursor_end_1based']), locus['strand'])
                row['equivalence_link'] = bool_str(locus['representative_precursor_sequence_DNA'].upper() == expected)
                row['parent_containment'] = bool_str(f['start'] <= region[1] and region[2] <= f['end'])
                notes.append(f'representative_model={locus["representative_model_id"]}')
                if (region[1], region[2]) != (f['start'], f['end']):
                    notes.append('gene_span_wider_than_representative_precursor')
            else:
                notes.append('locus_row_missing')
        elif f['type'] == 'miRNA_primary_transcript':
            parent = by_id.get(a.get('Parent', ''))
            row['parent_containment'] = bool_str(parent is not None and parent['start'] <= f['start'] and f['end'] <= parent['end'])
            locus = mir_loci.get(a.get('Parent', ''))
            if a.get('representative') == 'true' and locus:
                expected = prec.get(locus['release_locus_id'])
                row['expected_sequence_source'] = 'precursor FASTA of parent locus (representative model only)'
                row['expected_sequence_id'] = locus['release_locus_id']
                row['fasta_id'] = locus['release_locus_id']
                row['fasta_link'] = bool_str(expected is not None)
                coords_match = (int(locus['representative_precursor_start_1based']), int(locus['representative_precursor_end_1based'])) == (f['start'], f['end'])
                if not coords_match:
                    notes.append('representative_transcript_coords_differ_from_loci_table')
            else:
                row['expected_sequence_source'] = 'none (non-representative model; no per-model sequence in release)'
        elif f['type'] == 'miRNA':
            pid = a.get('reference_product_id', '')
            expected = mature.get(pid)
            prod = mir_products.get(pid)
            row['expected_sequence_source'] = 'miRNA_bona_fide.mature_and_star.fa via reference_product_id'
            row['expected_sequence_id'] = pid
            row['fasta_id'] = pid if expected else ''
            row['fasta_link'] = bool_str(expected is not None)
            parent = by_id.get(a.get('Parent', ''))
            row['parent_containment'] = bool_str(parent is not None and parent['start'] <= f['start'] and f['end'] <= parent['end'])
            if prod:
                locs = parse_locus_list(prod['genomic_locations'], genome)
                row['listed_placements'] = len(locs)
                row['placement_includes_feature'] = bool_str((f['seqid'], f['start'], f['end'], f['strand']) in locs)
                row['multimapper'] = bool_str(len(locs) > 1)
                row['equivalence_link'] = bool_str(prod['sequence_equivalence_group'] == a.get('sequence_equivalence_group')
                                                   and prod['mature_sequence_DNA'].upper() == expected)
                notes.append(f'role={a.get("product_role")};arm={a.get("arm")}')
            else:
                notes.append('product_row_missing')
        if expected is not None and row['coordinate_status'] == 'valid':
            fwd = genome.fetch(region[0], region[1], region[2])
            row['genome_match_orientation'] = orientation_match(expected, fwd)
            others = seqsets['miRNA_bona_fide.mature_and_star.fa'] if f['type'] == 'miRNA' else seqsets['miRNA_bona_fide.precursors.fa']
            row['sequence_status'] = sequence_status(row['genome_match_orientation'], 'transcript', region[3], fwd, others)
        else:
            row['sequence_status'] = 'not_evaluable'
        evaluable = row['sequence_status'] != 'not_evaluable'
        row['feature_pass'] = bool_str(
            row['coordinate_status'] == 'valid'
            and (row['sequence_status'] == 'exact_for_ID' if evaluable else True)
            and row['fasta_link'] in ('yes', '') and row['equivalence_link'] in ('yes', '')
            and row['parent_containment'] in ('yes', '') and row['placement_includes_feature'] in ('yes', '')
            and not any('differ' in n or 'missing' in n for n in notes if 'wider' not in n))
        row['notes'] = ';'.join(notes)
        rows.append(row)

    # ------------------------------------------------------------------ structural ncRNA
    exons_by_parent = collections.defaultdict(list)
    for f in feats:
        if f['type'] == 'exon' and f['attrs'].get('biotype') == 'tRNA_exon':
            exons_by_parent[f['attrs']['Parent']].append(f)
    struct_fasta = {'rRNA': 'rRNA.fa', 'snRNA': 'snRNA.fa', 'snoRNA': 'snoRNA.fa', 'tRNA': 'tRNA.genomic.fa'}
    single_pos_convention = collections.Counter()
    for f in feats:
        a = f['attrs']
        if f['type'] not in struct_fasta or a.get('biotype') != f['type']:
            continue
        cls = f['type']
        fk = struct_fasta[cls]
        row = new_row(f, cls, a.get('evidence_tier', ''))
        row['coordinate_status'] = coordinate_status(f, genome)
        row['expected_orientation'] = 'transcript (GFF strand)'
        notes = []
        member = member_by_id.get(a.get('source_feature_id', ''))
        expected = None
        row['expected_sequence_source'] = f'{fk} via structural_ncRNA.sequence_membership.tsv (member_id=source_feature_id)'
        if member:
            rid = member['reference_id']
            row['expected_sequence_id'] = rid
            expected = fasta[fk].get(rid)
            row['fasta_id'] = rid if expected else ''
            row['fasta_link'] = bool_str(expected is not None and expected == member['canonical_dna_sequence'].upper())
            eq = seqeq_by_id.get(rid)
            row['equivalence_link'] = bool_str(eq is not None and eq['sequence_DNA'].upper() == expected)
            locs = parse_locus_list(member['exact_B73_loci'], genome)
            row['listed_placements'] = len(locs)
            hit = (f['seqid'], f['start'], f['end'], f['strand']) in locs
            if not hit and any(e is None for (_, _, e, _) in locs):
                if (f['seqid'], f['start'], None, f['strand']) in locs:
                    hit = True
                    single_pos_convention['start_1based'] += 1
                elif (f['seqid'], f['start'] - 1, None, f['strand']) in locs:
                    hit = True
                    single_pos_convention['start_0based'] += 1
            row['placement_includes_feature'] = bool_str(hit)
            row['multimapper'] = bool_str(len(locs) > 1)
        else:
            row['fasta_link'] = 'no'
            notes.append('membership_row_missing')
        if expected is not None and row['coordinate_status'] == 'valid':
            if cls == 'tRNA':
                exons = sorted(exons_by_parent.get(a['ID'], []), key=lambda e: e['start'])
                notes.append(f'exons={len(exons)}')
                spliced = ''.join(genome.fetch(e['seqid'], e['start'], e['end']) for e in exons)
                full = genome.fetch(f['seqid'], f['start'], f['end'])
                m_spliced = orientation_match(expected, spliced) if exons else 'none'
                m_full = orientation_match(expected, full)
                if len(exons) > 1:
                    notes.append(f'spliced_match={m_spliced};full_span_match={m_full}')
                use = spliced if m_spliced != 'none' else full
                row['genome_match_orientation'] = m_spliced if m_spliced != 'none' else m_full
                if len(exons) > 1:
                    notes.append('sequence_model=' + ('intron_removed' if m_spliced != 'none' else
                                                      'full_span' if m_full != 'none' else 'neither'))
                fwd = use
            else:
                fwd = genome.fetch(f['seqid'], f['start'], f['end'])
                row['genome_match_orientation'] = orientation_match(expected, fwd)
            row['sequence_status'] = sequence_status(row['genome_match_orientation'], 'transcript', f['strand'],
                                                     fwd, seqsets[fk])
        else:
            row['sequence_status'] = 'not_evaluable'
        row['feature_pass'] = bool_str(row['coordinate_status'] == 'valid' and row['sequence_status'] == 'exact_for_ID'
                                       and row['fasta_link'] == 'yes' and row['equivalence_link'] == 'yes'
                                       and row['placement_includes_feature'] == 'yes')
        row['notes'] = ';'.join(notes)
        rows.append(row)
    for f in feats:
        if f['type'] == 'exon' and f['attrs'].get('biotype') == 'tRNA_exon':
            row = new_row(f, 'tRNA_exon', f['attrs'].get('evidence_tier', ''))
            row['coordinate_status'] = coordinate_status(f, genome)
            row['sequence_status'] = 'not_evaluable'
            row['expected_sequence_source'] = 'evaluated through parent tRNA'
            parent = by_id.get(f['attrs']['Parent'])
            row['parent_containment'] = bool_str(parent is not None and parent['start'] <= f['start'] and f['end'] <= parent['end'])
            row['feature_pass'] = bool_str(row['coordinate_status'] == 'valid' and row['parent_containment'] == 'yes')
            rows.append(row)

    # ------------------------------------------------------------------ exploratory MIR coordinate models
    for f in feats:
        a = f['attrs']
        if f['source'] != 'B73v5_MIR_candidate_review':
            continue
        row = new_row(f, 'unsupported_MIR_candidate_' + f['type'], a.get('evidence_tier', ''))
        row['coordinate_status'] = coordinate_status(f, genome)
        row['sequence_status'] = 'not_evaluable'
        row['expected_sequence_source'] = 'none (no per-model precursor sequence in release); products checked in exploratory_products.tsv'
        if f['type'] == 'ncRNA':
            parent = by_id.get(a.get('Parent', ''))
            row['parent_containment'] = bool_str(parent is not None and parent['start'] <= f['start'] and f['end'] <= parent['end'])
        row['feature_pass'] = bool_str(row['coordinate_status'] == 'valid' and row['parent_containment'] in ('yes', ''))
        rows.append(row)

    # sanity: every GFF feature got exactly one row
    covered = collections.Counter(r['feature_id'] for r in rows)
    unhandled = [f['attrs'].get('ID') for f in feats if covered[f['attrs'].get('ID')] != 1]
    if unhandled:
        notes_global.append(f'features_without_exactly_one_row={len(unhandled)} e.g. {unhandled[:3]}')

    # ------------------------------------------------------------------ exploratory products (FASTA-record level)
    gff_mir_models = collections.defaultdict(list)
    for f in feats:
        if f['source'] == 'B73v5_MIR_candidate_review' and f['type'] == 'ncRNA_gene':
            gff_mir_models[f['attrs'].get('source_internal_id', '')].append(f)
    phaslike_loci = {f['attrs']['Name']: f for f in feats
                     if f['type'] == 'ncRNA_gene' and f['attrs'].get('biotype', '').endswith('PHAS-like')}
    xrows = []
    for cls, fk, kind in (('unsupported_MIR_candidate', 'miRNA_unsupported_candidates.fa', 'mir'),
                          ('21-PHAS-like_candidate', '21_PHAS_like.candidate_products.fa', 'phas'),
                          ('24-PHAS-like_candidate', '24_PHAS_like.candidate_products.fa', 'phas')):
        for rid, s in fasta[fk].items():
            eq = seqeq_by_id.get(rid)
            members = eq['member_ids'].split(',') if eq else []
            regions = []
            for m in members:
                if kind == 'mir':
                    regions += gff_mir_models.get(m, [])
                elif m in phaslike_loci:
                    regions.append(phaslike_loci[m])
            fwd_hit = rc_hit = 0
            for g in regions:
                gs = genome.fetch(g['seqid'], g['start'], g['end'])
                strand_seq = rc(gs) if g['strand'] == '-' else gs
                if kind == 'mir':
                    if s in strand_seq:
                        fwd_hit += 1
                    elif rc(s) in strand_seq:
                        rc_hit += 1
                else:
                    if s in gs:
                        fwd_hit += 1
                    elif rc(s) in gs:
                        rc_hit += 1
            if not eq:
                status = 'unresolved_identifier'
            elif not regions:
                status = 'not_evaluable_no_coordinate_model'
            elif fwd_hit:
                status = 'contained_in_own_region_expected_orientation'
            elif rc_hit:
                status = 'contained_in_own_region_opposite_orientation'
            else:
                status = 'absent_from_own_regions'
            xrows.append({'class': cls, 'evidence_tier': 'exploratory', 'fasta_id': rid, 'sequence': s,
                          'member_ids': ','.join(members), 'regions_found': len(regions),
                          'regions_with_expected_orientation_hit': fwd_hit,
                          'regions_with_opposite_orientation_hit_only': rc_hit, 'status': status,
                          'expected_orientation': 'transcript (model strand)' if kind == 'mir' else 'reference_forward (PHAS-like loci are unstranded)'})

    # ------------------------------------------------------------------ FASTA coverage (canonical)
    used = collections.defaultdict(set)
    for r in rows:
        if r['fasta_id']:
            for fid in r['fasta_id'].split(','):
                used[r['class']].add(fid)
    cov_rows = []
    cov_spec = [('21-PHAS_phasiRNA', '21_PHAS.phasiRNAs.fa', ['21-PHAS_phasiRNA']),
                ('24-PHAS_phasiRNA', '24_PHAS.phasiRNAs.fa', ['24-PHAS_phasiRNA']),
                ('miRNA_mature_star', 'miRNA_bona_fide.mature_and_star.fa', ['miRNA_miRNA']),
                ('miRNA_precursor', 'miRNA_bona_fide.precursors.fa', ['miRNA_gene', 'miRNA_miRNA_primary_transcript']),
                ('tRNA', 'tRNA.genomic.fa', ['tRNA']), ('rRNA', 'rRNA.fa', ['rRNA']),
                ('snRNA', 'snRNA.fa', ['snRNA']), ('snoRNA', 'snoRNA.fa', ['snoRNA'])]
    for label, fk, classes in cov_spec:
        referenced = set().union(*[used[c] for c in classes])
        unref = sorted(set(fasta[fk]) - referenced)
        cov_rows.append({'fasta': fk, 'records': len(fasta[fk]), 'unique_sequences': len(seq_to_ids[fk]),
                         'records_referenced_by_GFF_features': len(set(fasta[fk]) & referenced),
                         'unreferenced_records': len(unref), 'unreferenced_examples': ','.join(unref[:5])})

    # CCA-aware consistency (not a genome check)
    cca = fasta['tRNA.CCA_aware_quantification.fa']
    gen = fasta['tRNA.genomic.fa']
    cca_stats = collections.Counter()
    for rid, s in gen.items():
        cid = rid.replace('_tRNA_', '_tRNA_CCA_')
        c = cca.get(cid)
        if c is None:
            cca_stats['missing_partner'] += 1
        elif s.endswith('CCA') and c == s:
            cca_stats['genomic_already_CCA_unchanged'] += 1
        elif not s.endswith('CCA') and c == s + 'CCA':
            cca_stats['CCA_appended'] += 1
        else:
            cca_stats['inconsistent'] += 1

    # ------------------------------------------------------------------ write outputs
    out = args.outdir

    def write_tsv(name, fields, data):
        with open(os.path.join(out, name), 'w', newline='') as handle:
            w = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(data)

    write_tsv('features.all_classes.tsv', ROW_FIELDS, rows)
    per_class = collections.defaultdict(list)
    for r in rows:
        per_class[r['class']].append(r)
    os.makedirs(os.path.join(out, 'per_class'), exist_ok=True)
    for cls, data in per_class.items():
        write_tsv(os.path.join('per_class', f'{cls}.features.tsv'), ROW_FIELDS, data)
    write_tsv('exploratory_products.tsv', list(xrows[0].keys()) if xrows else ['class'], xrows)
    write_tsv('fasta_coverage.tsv', list(cov_rows[0].keys()), cov_rows)

    summary = []
    for cls, data in sorted(per_class.items(), key=lambda kv: (kv[1][0]['evidence_tier'], kv[0])):
        st = collections.Counter(r['sequence_status'] for r in data)
        summary.append({
            'class': cls, 'evidence_tier': ','.join(sorted({r['evidence_tier'] for r in data})),
            'total_features': len(data),
            'compatible_coordinate': sum(r['coordinate_status'] == 'valid' for r in data),
            'sequence_compatible_any_orientation': st['exact_for_ID'] + st['reverse_complement_for_ID'],
            'exact_sequence_for_feature': st['exact_for_ID'],
            'reverse_complement_for_feature': st['reverse_complement_for_ID'],
            'only_other_ID_at_coordinate': st['only_other_ID_at_coordinate'],
            'no_match': st['absent'],
            'not_evaluable': st['not_evaluable'],
            'multimapper_features': sum(r['multimapper'] == 'yes' for r in data),
            'fasta_link_fail': sum(r['fasta_link'] == 'no' for r in data),
            'equivalence_link_fail': sum(r['equivalence_link'] == 'no' for r in data),
            'observation_link_fail': sum(r['observation_link'] not in ('yes', '') for r in data),
            'observation_counts_disagree': sum(r['observation_counts_agree'] == 'no' for r in data),
            'placement_not_listed': sum(r['placement_includes_feature'] == 'no' for r in data),
            'parent_containment_fail': sum(r['parent_containment'] == 'no' for r in data),
            'features_with_notes': sum(bool(r['notes']) for r in data),
            'feature_pass': sum(r['feature_pass'] == 'yes' for r in data),
            'feature_fail': sum(r['feature_pass'] == 'no' for r in data),
        })
    write_tsv('summary.per_class.tsv', list(summary[0].keys()), summary)
    xsum = collections.Counter((r['class'], r['status']) for r in xrows)
    write_tsv('summary.exploratory_products.tsv', ['class', 'status', 'count'],
              [{'class': c, 'status': s, 'count': n} for (c, s), n in sorted(xsum.items())])

    fails = [r for r in rows if r['feature_pass'] == 'no']
    examples = []
    seen = collections.Counter()
    for r in fails:
        if seen[r['class']] < 10:
            examples.append(r)
            seen[r['class']] += 1
    write_tsv('failing_examples.tsv', ROW_FIELDS, examples)

    # PHAS orientation cross-tab (strand x genome match)
    xtab = collections.Counter((r['class'], r['strand'], r['genome_match_orientation'])
                               for r in rows if r['class'].endswith('PHAS_phasiRNA'))
    write_tsv('phasiRNA_orientation_crosstab.tsv', ['class', 'gff_strand', 'genome_match_orientation', 'count'],
              [{'class': c, 'gff_strand': s, 'genome_match_orientation': m, 'count': n} for (c, s, m), n in sorted(xtab.items())])

    # ------------------------------------------------------------------ provenance
    sums = {}
    with open(gpath('provenance', 'SHA256SUMS')) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                sums[parts[-1].lstrip('*').lstrip('./')] = parts[0]
    checked = ['annotation/B73v5_sRNA_reference_v1.combined.gff3'] + [f'fasta/{k}' for k in fasta] + [
        'tables/sequence_equivalence.tsv', 'tables/structural_ncRNA.sequence_membership.tsv',
        'tables/miRNA_bona_fide.products.tsv', 'tables/miRNA_bona_fide.loci.tsv',
        'tables/21_PHAS.phasiRNA_observations.tsv.gz', 'tables/24_PHAS.phasiRNA_observations.tsv.gz']
    prov = []
    for rel in checked:
        h = sha256(gpath(rel))
        prov.append({'input': rel, 'sha256': h, 'release_SHA256SUMS': sums.get(rel, 'not_listed'),
                     'matches_release': bool_str(sums.get(rel) == h) if rel in sums else 'NA'})
    real_genome = os.path.realpath(args.genome)
    prov.append({'input': f'genome:{real_genome}',
                 'sha256': 'skipped' if args.skip_genome_sha else sha256(real_genome),
                 'release_SHA256SUMS': 'report: a40ee35888a9c0b956c448752f59637639a7906d7fa14971c9047eb088ded2c2',
                 'matches_release': 'see_report'})
    prov.append({'input': 'genome.fai', 'sha256': sha256(args.genome + '.fai'), 'release_SHA256SUMS': 'NA', 'matches_release': 'NA'})
    for phase, p in (('21', args.sam21), ('24', args.sam24)):
        if p:
            prov.append({'input': f'sam{phase}:{os.path.realpath(p)}', 'sha256': sha256(p),
                         'release_SHA256SUMS': 'NA (generated 2026-11-04)', 'matches_release': 'NA'})
    write_tsv('input_checksums.tsv', ['input', 'sha256', 'release_SHA256SUMS', 'matches_release'], prov)

    meta = {
        'script': os.path.realpath(__file__), 'script_sha256': sha256(os.path.realpath(__file__)),
        'command': ' '.join(sys.argv), 'run_utc': datetime.datetime.utcnow().isoformat(timespec='seconds') + 'Z',
        'python': platform.python_version(), 'host': platform.node(),
        'gff_features': len(feats), 'feature_rows': len(rows), 'exploratory_product_rows': len(xrows),
        'structural_single_position_convention': dict(single_pos_convention),
        'tRNA_CCA_aware_consistency': dict(cca_stats), 'global_notes': notes_global,
    }
    try:
        meta['samtools'] = subprocess.run(['samtools', '--version'], capture_output=True, text=True).stdout.splitlines()[0]
    except Exception:
        meta['samtools'] = 'unavailable'
    with open(os.path.join(out, 'run_metadata.json'), 'w') as handle:
        json.dump(meta, handle, indent=2)

    # console
    cols = ['class', 'evidence_tier', 'total_features', 'compatible_coordinate', 'exact_sequence_for_feature',
            'reverse_complement_for_feature', 'only_other_ID_at_coordinate', 'no_match', 'not_evaluable',
            'multimapper_features', 'feature_pass', 'feature_fail']
    print('\t'.join(cols))
    for s in summary:
        print('\t'.join(str(s[c]) for c in cols))
    print('\nexploratory products:')
    for (c, s), n in sorted(xsum.items()):
        print(f'  {c}\t{s}\t{n}')
    print('\nFASTA coverage:')
    for c in cov_rows:
        print('  ' + '\t'.join(str(v) for v in c.values()))
    print('\nmeta:', json.dumps({k: meta[k] for k in ('gff_features', 'feature_rows', 'structural_single_position_convention',
                                                     'tRNA_CCA_aware_consistency', 'global_notes')}))


if __name__ == '__main__':
    main()
