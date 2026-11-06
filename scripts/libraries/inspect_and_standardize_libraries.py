#!/usr/bin/env python3
"""Parse each downloaded GEO small-RNA file, report sanity statistics, and write a standard tag-count file.

Standard tag-count (Phasis -libformat T): SEQ<TAB>COUNT, no header, DNA alphabet (U->T), uppercase,
one line per distinct sequence, sorted by sequence. Only written when every count is a positive integer
(raw abundances). Sequences with characters other than ACGT, or outside --min-len..--max-len, are dropped
and counted. No adapter trimming or quality filtering is done (processed files are used as deposited).

Supported layouts (auto-detected, recorded): seq<TAB>count, count<TAB>seq (optional header line; extra
columns ignored after inspection), collapsed FASTA with count in the header (x123, (123), _123, count=123,
or header == sequence with trailing count), and FASTA of individual reads (each record counts 1).
"""
import argparse
import collections
import csv
import gzip
import hashlib
import os
import re
import sys

SEQ = re.compile(r'^[ACGTUN]+$', re.I)
NUM = re.compile(r'^\d+(\.\d+)?([eE][+-]?\d+)?$')
HEAD_COUNT = [re.compile(p) for p in (r'_x(\d+)$', r'x(\d+)$', r'\((\d+)\)$', r'count[=_:](\d+)', r'[_\-|:](\d+)$', r'\s(\d+)$')]
ADAPTERS = {'illumina_3p_TGGAATTCTCGG': 'TGGAATTCTCGG', 'illumina_truseq_AGATCGGAAGAGC': 'AGATCGGAAGAGC',
            'nextflex_3p_CTGTAGGCAC': 'CTGTAGGCAC', 'polyA_AAAAAAAAAA': 'AAAAAAAAAA'}


def opener(path):
    return gzip.open(path, 'rt', errors='replace') if path.endswith('.gz') else open(path, errors='replace')


def parse(path):
    """Return (counts dict raw_seq->value, layout, notes)."""
    counts = collections.defaultdict(float)
    notes = []
    with opener(path) as fh:
        first = ''
        for first in fh:
            if first.strip():
                break
        if first.startswith('>'):
            header, seq_parts, layout = first[1:].strip(), [], None
            n_counted = n_records = 0

            def flush(h, s):
                nonlocal layout, n_counted, n_records
                if not s:
                    return
                n_records += 1
                c = None
                for p in HEAD_COUNT:
                    m = p.search(h)
                    if m:
                        c = float(m.group(1)); break
                if c is None:
                    c = 1.0
                else:
                    n_counted += 1
                counts[s.upper()] += c
            for line in fh:
                if line.startswith('>'):
                    flush(header, ''.join(seq_parts)); header, seq_parts = line[1:].strip(), []
                else:
                    seq_parts.append(line.strip())
            flush(header, ''.join(seq_parts))
            layout = 'fasta_collapsed' if n_counted >= 0.99 * n_records else ('fasta_reads' if n_counted == 0 else 'fasta_mixed')
            if layout == 'fasta_mixed':
                notes.append(f'{n_counted}/{n_records} headers carried counts')
            return counts, layout, notes
        lines = [first] + [next(fh, '') for _ in range(5)]
        fields = [re.split(r'[\t,]+', l.strip()) for l in lines if l.strip()]
        header = None
        if fields and not (SEQ.match(fields[0][0]) or NUM.match(fields[0][0])):
            header = lines[0].strip(); fields = fields[1:]
        if not fields:
            return counts, 'empty', notes
        if SEQ.match(fields[0][0]) and len(fields[0]) > 1 and NUM.match(fields[0][1]):
            si, ci, layout = 0, 1, 'seq_count'
        elif NUM.match(fields[0][0]) and len(fields[0]) > 1 and SEQ.match(fields[0][1]):
            si, ci, layout = 1, 0, 'count_seq'
        elif SEQ.match(fields[0][0]) and len(fields[0]) == 1:
            si, ci, layout = 0, None, 'seq_only_lines'
        else:
            return counts, 'unparsed', [f'first fields: {fields[0][:4]}']
        if header:
            notes.append('header=' + header[:80].replace('\t', '|'))
        ncols = collections.Counter()
        for line in lines[(1 if header else 0):] + list(fh):
            f = re.split(r'[\t,]+', line.strip())
            if not f or not f[0]:
                continue
            ncols[len(f)] += 1
            try:
                s = f[si]
                c = 1.0 if ci is None else float(f[ci])
            except (IndexError, ValueError):
                notes.append('unparsable_line'); continue
            counts[s.upper()] += c
        if max(ncols, default=0) > 2:
            notes.append(f'extra_columns_ignored:{dict(ncols)}')
        return counts, layout, notes


def process(job):
    m, a, pool = job
    path = os.path.join(os.path.dirname(a.manifest), m['gse'], m['url'].split('/')[-1])
    raw, layout, notes = parse(path)
    st = collections.Counter()
    lens = collections.Counter()
    five = {21: collections.Counter(), 24: collections.Counter()}
    adapt = collections.Counter()
    has_u = any('U' in s for s in list(raw)[:20000])
    integer = all(float(v).is_integer() and v > 0 for v in raw.values())
    clean = collections.Counter()
    for s, c in raw.items():
        st['raw_reads'] += c
        st['raw_distinct'] += 1
        d = s.replace('U', 'T')
        if not re.fullmatch(r'[ACGT]+', d):
            st['dropped_nonACGT_reads'] += c; continue
        if not (a.min_len <= len(d) <= a.max_len):
            st['dropped_length_reads'] += c; continue
        clean[d] += c
        lens[len(d)] += c
        if len(d) in five:
            five[len(d)][d[0]] += c
        for k, ad in ADAPTERS.items():
            if ad in d:
                adapt[k] += c
    tot = sum(clean.values()) or 1
    out = os.path.join(a.outdir, f"{m['gsm']}.tag")
    written = ''
    if integer and clean:
        with open(out, 'w') as fh:
            for s in sorted(clean):
                fh.write(f'{s}\t{int(clean[s])}\n')
        written = hashlib.sha256(open(out, 'rb').read()).hexdigest()
    p = pool.get(m['gsm'], {})
    row = {'gsm': m['gsm'], 'gse': m['gse'], 'title': p.get('title', ''), 'organ': p.get('organ', ''),
           'subspecies': p.get('subspecies', ''), 'source_file': os.path.basename(path), 'source_sha256': m['sha256'],
           'layout': layout, 'integer_counts': integer, 'has_U': has_u,
           'raw_reads': int(st['raw_reads']), 'raw_distinct': st['raw_distinct'],
           'kept_reads': int(sum(clean.values())), 'kept_distinct': len(clean),
           'dropped_nonACGT_reads': int(st['dropped_nonACGT_reads']), 'dropped_length_reads': int(st['dropped_length_reads']),
           'min_len': min(lens) if lens else '', 'max_len': max(lens) if lens else '',
           'frac_18_30': round(sum(v for l, v in lens.items() if 18 <= l <= 30) / tot, 4),
           'frac_21': round(lens[21] / tot, 4), 'frac_22': round(lens[22] / tot, 4), 'frac_24': round(lens[24] / tot, 4),
           'top_length': lens.most_common(1)[0][0] if lens else '',
           'five_prime_21': ','.join(f'{n}:{round(five[21][n] / max(1, sum(five[21].values())), 2)}' for n in 'ACGT'),
           'five_prime_24': ','.join(f'{n}:{round(five[24][n] / max(1, sum(five[24].values())), 2)}' for n in 'ACGT'),
           'adapter_traces': ';'.join(f'{k}:{round(v / tot, 4)}' for k, v in adapt.items() if v / tot >= 0.001),
           'tagcount_file': os.path.basename(out) if written else '', 'tagcount_sha256': written,
           'notes': ';'.join(notes)}
    print(m['gsm'], layout, row['kept_reads'], row['integer_counts'], file=sys.stderr, flush=True)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--pool', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--stats', required=True)
    ap.add_argument('--min-len', type=int, default=15)
    ap.add_argument('--max-len', type=int, default=40)
    ap.add_argument('--jobs', type=int, default=1)
    a = ap.parse_args()
    pool = {r['gsm']: r for r in csv.DictReader(open(a.pool), delimiter='\t')}
    os.makedirs(a.outdir, exist_ok=True)
    jobs = [(m, a, pool) for m in csv.DictReader(open(a.manifest), delimiter='\t') if m['status'] == 'ok']
    import multiprocessing
    with multiprocessing.Pool(a.jobs) as mp:
        rows = mp.map(process, jobs, chunksize=1)
    with open(a.stats, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(rows)


if __name__ == '__main__':
    main()
