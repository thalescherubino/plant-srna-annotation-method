#!/usr/bin/env python3
"""Split / merge miRScore runs for parallel execution (2026-11-07, user request: 8 chunks x 4 cores).

split: <set_dir>/inputs/hairpin.fa is dealt round-robin into N chunk hairpin files
       (<set_dir>/chunks/cNN/inputs/hairpin.fa). Every chunk gets the full mature.fa: miRScore only looks mature
       records up by hairpin name, so unused mature records are ignored. Each hairpin is scored independently and
       bowtie -a reports every alignment of a read to every hairpin, so per-hairpin results do not depend on which
       other hairpins share the run.
merge: concatenates each chunk's <name>_miRScore_results.csv, _reads.csv, _alt_mirna_results.csv and _alt_reads.csv
       into <set_dir>/merged/<set>_*.csv (header once), checks every input hairpin appears exactly once in the
       merged results, and writes merged/MERGE_SUMMARY.txt.
"""
import argparse
import csv
import glob
import os
import sys


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def split(set_dir, n):
    hp = list(read_fasta(os.path.join(set_dir, 'inputs', 'hairpin.fa')))
    mature = os.path.join(set_dir, 'inputs', 'mature.fa')
    for i in range(n):
        d = os.path.join(set_dir, 'chunks', f'c{i + 1:02d}', 'inputs')
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'hairpin.fa'), 'w') as out:
            out.writelines(f'>{h}\n{s}\n' for h, s in hp[i::n])
        dst = os.path.join(d, 'mature.fa')
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(os.path.abspath(mature), dst)
    print(f'{set_dir}: {len(hp)} hairpins -> {n} chunks of {len(hp) // n}-{-(-len(hp) // n)}')


def merge(set_dir, name):
    chunks = sorted(glob.glob(os.path.join(set_dir, 'chunks', 'c*')))
    out_dir = os.path.join(set_dir, 'merged')
    os.makedirs(out_dir, exist_ok=True)
    summary = []
    for suffix in ('miRScore_results', 'reads', 'alt_mirna_results', 'alt_reads'):
        rows, header = [], None
        for c in chunks:
            f = os.path.join(c, 'results', f'{name}_{os.path.basename(c)}_{suffix}.csv')
            if not os.path.exists(f):
                summary.append(f'missing {f}')
                continue
            with open(f, newline='') as fh:
                r = list(csv.reader(fh))
            if not r:
                continue
            if header is None:
                header = r[0]
            elif r[0] != header:
                sys.exit(f'header mismatch in {f}')
            rows += r[1:]
        if header is not None:
            with open(os.path.join(out_dir, f'{name}_{suffix}.csv'), 'w', newline='') as out:
                w = csv.writer(out)
                w.writerow(header)
                w.writerows(rows)
        summary.append(f'{suffix}: {len(rows)} rows from {len(chunks)} chunks')
    want = [h for h, _ in read_fasta(os.path.join(set_dir, 'inputs', 'hairpin.fa'))]
    with open(os.path.join(out_dir, f'{name}_miRScore_results.csv'), newline='') as fh:
        got = [r['name'] for r in csv.DictReader(fh)]
    missing = sorted(set(want) - set(got))
    dup = sorted({g for g in got if got.count(g) > 1})
    summary.append(f'hairpins in input: {len(want)}; in merged results: {len(got)}; missing: {len(missing)}; duplicated: {len(dup)}')
    if missing:
        summary.append('missing ids: ' + ' '.join(missing[:50]))
    with open(os.path.join(out_dir, 'MERGE_SUMMARY.txt'), 'w') as out:
        out.writelines(s + '\n' for s in summary)
    print('\n'.join(summary))
    if missing or dup:
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['split', 'merge'])
    ap.add_argument('set_dir')
    ap.add_argument('--chunks', type=int, default=8)
    ap.add_argument('--name', help='set name used for miRScore -n (merge)')
    a = ap.parse_args()
    if a.mode == 'split':
        split(a.set_dir, a.chunks)
    else:
        merge(a.set_dir, a.name or os.path.basename(os.path.normpath(a.set_dir)))


if __name__ == '__main__':
    main()
