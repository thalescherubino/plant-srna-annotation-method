#!/usr/bin/env python3
"""Per-library QC table and aggregate summary for a library selection built from raw SRA (added 2026-11-09, tomato).

Inputs: the selection prefix written by select_libraries.py (<prefix>.selected.tsv, <prefix>.excluded.tsv), the
library_stats.tsv and the <SRX>.report.json files written by sra_to_tagcount.py.
Adds review-only columns that are not hard filters:
  top_seq_share   fraction of kept reads from the single most abundant sequence (a library dominated by one
                  structural-RNA fragment or one miRNA shows a skewed size profile and 5' composition);
  soft_flags      21nt_5U_low (5'U < 0.40 among 21-nt reads), 24nt_5A_low (5'A < 0.40 among 24-nt reads),
                  dominated_by_one_seq (top_seq_share >= 0.20), trimmed_extra_bases (any extra 5'/3' cut).
Writes <out>.tsv (one row per processed library: selected or excluded, with the reason) and prints a summary.

Usage: library_qc_summary.py --prefix inventory/sra/selection_v1 --libdir libraries/sra --out inventory/sra/selection_v1.qc
"""
import argparse
import collections
import csv
import json
import os


def five(s, base):
    return float(dict(x.split(':') for x in s.split(','))[base])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--libdir', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    stats = {s['gsm']: s for s in csv.DictReader(open(os.path.join(a.libdir, 'library_stats.tsv')), delimiter='\t')}
    sel = {r['gsm']: r for r in csv.DictReader(open(a.prefix + '.selected.tsv'), delimiter='\t')}
    exc = {r['gsm']: r for r in csv.DictReader(open(a.prefix + '.excluded.tsv'), delimiter='\t') if r['gsm'] in stats}
    rows = []
    for g, s in stats.items():
        rep = json.load(open(os.path.join(a.libdir, 'tag', g + '.report.json')))
        with open(s['tagcount_file']) as fh:
            top = int(fh.readline().split('\t')[1])
        r = sel.get(g) or exc.get(g) or {}
        ch = rep['chosen']
        u21, a24 = five(s['five_prime_21'], 'T'), five(s['five_prime_24'], 'A')
        share = top / max(1, int(s['kept_reads']))
        soft = [f for f, c in (('21nt_5U_low', u21 < 0.4), ('24nt_5A_low', a24 < 0.4),
                               ('dominated_by_one_seq', share >= 0.2),
                               ('trimmed_extra_bases', ch['cut5'] or ch['cut3'])) if c]
        h = rep['length_hist']
        tot = max(1, sum(h.values()))
        rows.append({
            'library': g, 'status': 'selected' if g in sel else 'excluded', 'reason': '' if g in sel else r.get('exclusion_reason', 'not_selected'),
            'series': r.get('gse', ''), 'organ': r.get('organ', ''), 'genotype': r.get('subspecies', ''),
            'treatment': r.get('treatment_class', ''), 'adapter': ch['adapter'], 'cut5': ch['cut5'], 'cut3': ch['cut3'],
            'raw_M': round(int(s['raw_reads']) / 1e6, 2), 'kept_M': round(int(s['kept_reads']) / 1e6, 2),
            'keep_frac': round(int(s['kept_reads']) / max(1, int(s['raw_reads'])), 3),
            'main_peak': max(h, key=lambda k: h[k]), 'pct21': round(100 * h['21'] / tot, 1), 'pct22': round(100 * h['22'] / tot, 1),
            'pct23': round(100 * h['23'] / tot, 1), 'pct24': round(100 * h['24'] / tot, 1),
            '5U_21nt': u21, '5A_24nt': a24, 'exact_map_sample': float(s['map_rate_total_sample']),
            'top_seq_share': round(share, 3), 'qc_flags': s['flags'], 'soft_flags': ';'.join(soft)})
    rows.sort(key=lambda x: (x['status'] != 'selected', x['series'], x['library']))
    with open(a.out + '.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    S = [x for x in rows if x['status'] == 'selected']
    E = [x for x in rows if x['status'] == 'excluded']
    med = lambda v: sorted(v)[len(v) // 2] if v else 0  # noqa: E731
    print(f'processed={len(rows)} selected={len(S)} excluded_after_processing={len(E)}')
    print(f"selected reads: {sum(x['kept_M'] for x in S):.0f} M kept; median {med([x['kept_M'] for x in S])} M per library")
    print('selected by organ:', dict(collections.Counter(x['organ'] for x in S).most_common()))
    print('selected by genotype:', dict(collections.Counter(x['genotype'] for x in S).most_common()))
    print('selected by treatment:', dict(collections.Counter(x['treatment'] for x in S).most_common()))
    print('selected series:', len({x['series'] for x in S}), dict(collections.Counter(x['series'] for x in S).most_common()))
    print('selected by trimming:', dict(collections.Counter(f"{x['adapter']}/cut5={x['cut5']}/cut3={x['cut3']}" for x in S).most_common()))
    print('selected main peak:', dict(collections.Counter(x['main_peak'] for x in S).most_common()))
    print(f"selected medians: 21nt {med([x['pct21'] for x in S])}%, 24nt {med([x['pct24'] for x in S])}%, "
          f"5'U@21 {med([x['5U_21nt'] for x in S])}, 5'A@24 {med([x['5A_24nt'] for x in S])}, "
          f"exact map {med([x['exact_map_sample'] for x in S])}")
    print('selected soft flags:', dict(collections.Counter(f for x in S for f in x['soft_flags'].split(';') if f).most_common()))
    print('excluded after processing, by reason:', dict(collections.Counter(x['reason'] for x in E).most_common()))


if __name__ == '__main__':
    main()
