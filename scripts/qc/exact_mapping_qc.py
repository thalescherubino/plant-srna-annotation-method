#!/usr/bin/env python3
"""Exact-match genome mapping rate per library (sanity check; not an annotation step).

For each standardized tag-count library, distinct 18-30-nt sequences are aligned with Bowtie 1
(-v 0 -k 1, both strands) against the genome; read-weighted mapping fractions are reported overall and
for 21, 22 and 24 nt. Low rates for non-reference cultivars quantify reads lost to exact-match policies
(Phasis -mismat 0) because of genotype differences.
"""
import argparse
import collections
import csv
import os
import subprocess
import sys
import tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--index', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', type=int, default=32)
    ap.add_argument('--bowtie', default='bowtie')
    ap.add_argument('--tmpdir', default=None, help='scratch directory for temporary FASTA files')
    a = ap.parse_args()
    libdir = os.path.dirname(a.manifest)
    done = set()
    if os.path.exists(a.out):
        done = {r['library_id'] for r in csv.DictReader(open(a.out), delimiter='\t')}
    fields = ['library_id', 'gse', 'subspecies', 'organ', 'reads_18_30', 'mapped_18_30', 'frac_mapped_18_30',
              'frac_mapped_21', 'frac_mapped_22', 'frac_mapped_24', 'distinct_18_30', 'frac_distinct_mapped'] + \
             [f'reads_{n}' for n in range(18, 31)] + [f'frac_mapped_{n}nt' for n in range(18, 31)]
    fh = open(a.out, 'a', newline='')
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    if not done:
        w.writeheader()
    for m in csv.DictReader(open(a.manifest), delimiter='\t'):
        if m['library_id'] in done:
            continue
        counts, lens = [], []
        with tempfile.TemporaryDirectory(dir=a.tmpdir) as td:
            fa = os.path.join(td, 'q.fa')
            with open(fa, 'w') as out, open(os.path.join(libdir, m['tagcount_file'])) as tag:
                for i, line in enumerate(tag):
                    s, c = line.rstrip('\n').split('\t')
                    if 18 <= len(s) <= 30:
                        out.write(f'>{len(counts)}\n{s}\n')
                        counts.append(int(c)); lens.append(len(s))
            proc = subprocess.run([a.bowtie, '-f', '-v', '0', '-k', '1', '-p', str(a.threads), '--suppress', '2,3,4,5,6,7,8',
                                   a.index, fa], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            mapped = {int(x) for x in proc.stdout.split()}
        tot = collections.Counter(); mp = collections.Counter()
        for i, (c, l) in enumerate(zip(counts, lens)):
            tot['all'] += c; tot[l] += c
            if i in mapped:
                mp['all'] += c; mp[l] += c
        f = lambda k: round(mp[k] / tot[k], 4) if tot[k] else ''
        w.writerow({'library_id': m['library_id'], 'gse': m['gse'], 'subspecies': m['subspecies'], 'organ': m['organ'],
                    'reads_18_30': tot['all'], 'mapped_18_30': mp['all'], 'frac_mapped_18_30': f('all'),
                    'frac_mapped_21': f(21), 'frac_mapped_22': f(22), 'frac_mapped_24': f(24),
                    'distinct_18_30': len(counts), 'frac_distinct_mapped': round(len(mapped) / max(1, len(counts)), 4),
                    **{f'reads_{n}': tot[n] for n in range(18, 31)},
                    **{f'frac_mapped_{n}nt': f(n) for n in range(18, 31)}})
        fh.flush()
        print(m['library_id'], f('all'), file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
