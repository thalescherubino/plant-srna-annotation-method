#!/usr/bin/env python3
"""Build a spliced, CCA-aware mature tRNA quantification FASTA as a patch on top of
the frozen B73v5 curated sRNA reference v1 (2026-08-27).

NEW patch code (2026-11-05); not original release-building code. Read-only on the
release. The v1 GFF3 and tRNA.genomic.fa are NOT changed: they remain the
coordinate-bearing genomic (pre-tRNA, intron-containing) representation.

Problem fixed: for intron-containing tRNAs, v1 tRNA.genomic.fa holds the full genomic
span (5' exon + intron + 3' exon), and tRNA.CCA_aware_quantification.fa inherited that,
so mature (spliced) tRNA reads crossing the exon-exon junction could not map exactly.

Rules:
  * mature sequence = exons concatenated in transcript order (intron removed), taken
    from the genome at the v1 GFF3 exon coordinates;
  * append one terminal CCA only when the spliced sequence does not already end in CCA
    (the documented v1 rule);
  * identical mature sequences are collapsed into one record (sequence-equivalence);
  * a tRNA is spliced only if its intron call is canonical (anticodon ends exactly 1 nt
    before the intron, i.e. intron between positions 37/38) and its tRNAscan-SE score is
    >= --min-score; otherwise its v1 record is kept unchanged and flagged
    held_for_review (override with --splice-all).
"""
import argparse
import collections
import csv
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'validation_2026-11-05'))
import validate_srna_reference_strict as v  # noqa: E402  (genome reader, GFF/FASTA parsers)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def wrap(seq, width=80):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--release', required=True)
    ap.add_argument('--genome', required=True, help='faidx-indexed B73v5 FASTA')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--min-score', type=float, default=50.0,
                    help='minimum tRNAscan-SE score to splice automatically (default 50)')
    ap.add_argument('--splice-all', action='store_true', help='splice every multi-exon tRNA regardless of flags')
    args = ap.parse_args()

    R = args.release
    g = v.IndexedFasta(args.genome)
    gff_path = os.path.join(R, 'annotation', 'B73v5_sRNA_reference_v1.combined.gff3')
    feats = v.read_gff(gff_path)
    gen = v.read_fasta(os.path.join(R, 'fasta', 'tRNA.genomic.fa'))
    cca_path = os.path.join(R, 'fasta', 'tRNA.CCA_aware_quantification.fa')
    cca = v.read_fasta(cca_path)
    member = {r['member_id']: r for r in v.read_tsv(os.path.join(R, 'tables', 'structural_ncRNA.sequence_membership.tsv'))}

    # raw GFF score column is not kept by read_gff; re-read it for tRNA lines
    scores = {}
    with open(gff_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) == 9 and f[2] == 'tRNA':
                scores[v.parse_attrs(f[8])['ID']] = float(f[5]) if f[5] != '.' else float('nan')

    exons = collections.defaultdict(list)
    for f in feats:
        if f['type'] == 'exon' and f['attrs'].get('biotype') == 'tRNA_exon':
            exons[f['attrs']['Parent']].append(f)

    per_ref = collections.defaultdict(list)  # reference_id -> per-feature records
    for f in feats:
        if f['type'] != 'tRNA':
            continue
        a = f['attrs']
        rid = member[a['source_feature_id']]['reference_id']
        es = sorted(exons[a['ID']], key=lambda e: e['start'])
        full = g.fetch(f['seqid'], f['start'], f['end'])
        spliced = ''.join(g.fetch(e['seqid'], e['start'], e['end']) for e in es)
        if f['strand'] == '-':
            full, spliced = v.rc(full), v.rc(spliced)
        if full != gen[rid]:
            raise SystemExit(f'v1 genomic sequence mismatch for {a["ID"]} / {rid}')
        rec = {'feature': f, 'rid': rid, 'n_exons': len(es), 'full': full, 'spliced': spliced,
               'score': scores[a['ID']], 'introns': [], 'exon1_len': None, 'canonical': None}
        if len(es) > 1:
            first = es[0] if f['strand'] == '+' else es[-1]
            L = first['end'] - first['start'] + 1
            rec['exon1_len'] = L
            rec['canonical'] = spliced[L - 4:L - 1] == a.get('anticodon')
            rec['introns'] = [(es[i]['end'] + 1, es[i + 1]['start'] - 1) for i in range(len(es) - 1)]
        per_ref[rid].append(rec)

    # decide per reference_id (all members of a reference share sequence; verify structure agrees)
    decisions = {}
    for rid, recs in per_ref.items():
        multi = [r for r in recs if r['n_exons'] > 1]
        if not multi:
            continue
        if len(multi) != len(recs) or len({r['spliced'] for r in recs}) != 1:
            raise SystemExit(f'{rid}: members disagree on exon structure; refusing to patch')
        reasons = []
        if not all(r['canonical'] for r in recs):
            reasons.append('noncanonical_intron_position')
        low = min(r['score'] for r in recs)
        if low < args.min_score:
            reasons.append(f'low_tRNAscan_score({low:g}<{args.min_score:g})')
        decisions[rid] = ('spliced' if (args.splice_all or not reasons) else 'held_for_review', ';'.join(reasons))

    # build the new mature quantification set
    mature_of = {}  # genomic reference_id -> mature sequence (with CCA rule)
    for rid, gseq in gen.items():
        if decisions.get(rid, ('unchanged',))[0] == 'spliced':
            s = per_ref[rid][0]['spliced']
            mature_of[rid] = s if s.endswith('CCA') else s + 'CCA'
        else:
            old = cca[rid.replace('_tRNA_', '_tRNA_CCA_')]
            mature_of[rid] = old
    groups = collections.OrderedDict()
    for rid in gen:  # preserve v1 order
        groups.setdefault(mature_of[rid], []).append(rid)

    os.makedirs(os.path.join(args.outdir, 'fasta'), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, 'tables'), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, 'annotation'), exist_ok=True)

    fasta_out = os.path.join(args.outdir, 'fasta', 'tRNA.CCA_aware_quantification.spliced.fa')
    eq_rows, spl_n = [], 0
    with open(fasta_out, 'w') as out:
        for seq, rids in groups.items():
            touched = [r for r in rids if decisions.get(r, ('',))[0] == 'spliced']
            if touched:
                spl_n += 1
                qid = f'zma_B73v5_tRNA_CCA_SPL_{spl_n:04d}'
                status = 'spliced_mature_CCA_aware'
            else:
                assert len(rids) == 1, rids  # v1 genomic sequences are unique
                qid = rids[0].replace('_tRNA_', '_tRNA_CCA_')
                status = 'held_for_review_v1_unchanged' if rids[0] in decisions else 'v1_unchanged'
            feats_in = [r['feature']['attrs']['ID'] for rid in rids for r in per_ref[rid]]
            iso = sorted({r['feature']['attrs'].get('isotype', '') + r['feature']['attrs'].get('anticodon', '')
                          for rid in rids for r in per_ref[rid]})
            out.write(f'>{qid} status={status} genomic_refs={len(rids)} genomic_features={len(feats_in)} '
                      f'isotype_anticodon={",".join(iso)}\n{wrap(seq)}\n')
            eq_rows.append({'quantification_id': qid, 'status': status, 'mature_sequence_DNA': seq,
                            'length_nt': len(seq), 'ends_with_appended_or_encoded_CCA': seq.endswith('CCA'),
                            'genomic_reference_ids': ','.join(rids), 'genomic_feature_count': len(feats_in),
                            'genomic_feature_ids': ','.join(feats_in), 'isotype_anticodon': ','.join(iso)})

    def write_tsv(path, rows):
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)

    write_tsv(os.path.join(args.outdir, 'tables', 'tRNA_mature_sequence_equivalence.tsv'), eq_rows)

    patch_rows, intron_lines = [], []
    qid_of = {rid: row['quantification_id'] for row in eq_rows for rid in row['genomic_reference_ids'].split(',')}
    for rid, (decision, reason) in sorted(decisions.items()):
        for r in per_ref[rid]:
            f, a = r['feature'], r['feature']['attrs']
            patch_rows.append({
                'tRNA_feature_id': a['ID'], 'name': a.get('Name'), 'isotype': a.get('isotype'),
                'anticodon': a.get('anticodon'), 'seqid': f['seqid'], 'start_1based': f['start'],
                'end_1based': f['end'], 'strand': f['strand'], 'tRNAscan_score': r['score'],
                'introns_1based': ','.join(f'{s}-{e}' for s, e in r['introns']),
                'intron_length_nt': sum(e - s + 1 for s, e in r['introns']),
                'exon1_length_transcript_orientation': r['exon1_len'],
                'intron_canonical_after_anticodon': r['canonical'],
                'genomic_reference_id': rid, 'genomic_length_nt': len(r['full']),
                'spliced_length_nt': len(r['spliced']), 'decision': decision, 'flag_reasons': reason,
                'quantification_id': qid_of[rid],
                'genomic_sequence': r['full'], 'spliced_sequence': r['spliced'],
            })
            for i, (s, e) in enumerate(r['introns'], 1):
                attrs = (f'ID={a["ID"]}.intron{i};Parent={a["ID"]};biotype=tRNA_intron;evidence_tier=canonical;'
                         f'derived_from=v1_exon_gaps;mature_quantification_status={decision}')
                intron_lines.append('\t'.join([f['seqid'], 'tRNAscan-SE_2.0.13', 'intron', str(s), str(e),
                                               '.', f['strand'], '.', attrs]))
    write_tsv(os.path.join(args.outdir, 'tables', 'tRNA_intron_patch.tsv'), patch_rows)
    with open(os.path.join(args.outdir, 'annotation', 'B73v5_sRNA_reference_v1.tRNA_introns.sidecar.gff3'), 'w') as fh:
        fh.write('##gff-version 3\n##assembly Zm-B73-REFERENCE-NAM-5.0\n'
                 '# Sidecar to B73v5_sRNA_reference_v1.combined.gff3 (unchanged). Explicit introns of '
                 'multi-exon tRNAs, derived from gaps between v1 tRNA exon children.\n')
        fh.write('\n'.join(intron_lines) + '\n')

    # self-checks
    out_recs = v.read_fasta(fasta_out)
    assert len(set(out_recs.values())) == len(out_recs), 'duplicate sequences in output'
    for row in eq_rows:
        for rid in row['genomic_reference_ids'].split(','):
            if decisions.get(rid, ('',))[0] == 'spliced':
                core = per_ref[rid][0]['spliced']
                assert row['mature_sequence_DNA'] in (core, core + 'CCA')
            else:
                assert row['mature_sequence_DNA'] == cca[rid.replace('_tRNA_', '_tRNA_CCA_')]
    dec = collections.Counter(d for d, _ in decisions.values())
    feats_dec = collections.Counter(d for rid, (d, _) in decisions.items() for _ in per_ref[rid])
    summary = [
        ('v1_genomic_reference_sequences', len(gen)),
        ('v1_CCA_aware_records', len(cca)),
        ('multi_exon_tRNA_features', sum(feats_dec.values())),
        ('multi_exon_reference_sequences', len(decisions)),
        ('reference_sequences_spliced', dec['spliced']),
        ('features_spliced', feats_dec['spliced']),
        ('reference_sequences_held_for_review', dec['held_for_review']),
        ('features_held_for_review', feats_dec['held_for_review']),
        ('distinct_spliced_mature_sequences', spl_n),
        ('patched_CCA_aware_records', len(out_recs)),
        ('min_score_threshold', args.min_score),
        ('splice_all', args.splice_all),
        ('input_sha256:combined.gff3', sha256(gff_path)),
        ('input_sha256:tRNA.genomic.fa', sha256(os.path.join(R, 'fasta', 'tRNA.genomic.fa'))),
        ('input_sha256:tRNA.CCA_aware_quantification.fa', sha256(cca_path)),
        ('output_sha256:tRNA.CCA_aware_quantification.spliced.fa', sha256(fasta_out)),
        ('command', ' '.join(sys.argv)),
    ]
    with open(os.path.join(args.outdir, 'tables', 'patch_summary.tsv'), 'w') as fh:
        fh.write('metric\tvalue\n' + ''.join(f'{k}\t{val}\n' for k, val in summary))
    for k, val in summary:
        print(f'{k}\t{val}')


if __name__ == '__main__':
    main()
