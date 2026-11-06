#!/usr/bin/env python3
"""Probe the format of GEO per-sample supplementary files without full downloads.

For each candidate sample, stream the first bytes of each supplementary file
(HTTP range request), decompress what is available, and classify the content:

  tagcount_seq_count   SEQ<TAB>COUNT (Phasis -libformat T layout), optional header
  tagcount_count_seq   COUNT<TAB>SEQ
  fasta_collapsed      FASTA whose headers carry a read count
  fasta_reads          FASTA without counts (one record per read or unique sequence)
  fastq                FASTQ
  table_named          rows keyed by miRNA/gene/locus names (not read-level)
  sam_or_alignment     SAM/BAM/alignment-like text
  unknown

Read-only; writes inventory/<dir>/format_probe.tsv with the first lines kept for review.
"""
import argparse
import csv
import gzip
import io
import re
import sys
import time
import urllib.request
import zlib

UA = {'User-Agent': 'plant-srna-annotation-inventory/0.1'}
SEQ = re.compile(r'^[ACGTUNacgtun]{12,60}$')


def head_bytes(url, n=65536):
    url = url.replace('ftp://ftp.ncbi.nlm.nih.gov', 'https://ftp.ncbi.nlm.nih.gov')
    req = urllib.request.Request(url, headers={**UA, 'Range': f'bytes=0-{n - 1}'})
    for k in range(4):
        try:
            time.sleep(0.35)
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception:  # noqa: BLE001
            if k == 3:
                raise
            time.sleep(3 * (k + 1))


def decompress_partial(raw, name):
    if name.endswith('.gz'):
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            return d.decompress(raw)
        except zlib.error:
            return b''
    if name.endswith('.bz2'):
        return b''
    return raw


def classify(lines):
    lines = [l for l in lines if l.strip()][:200]
    if not lines:
        return 'unknown', ''
    if lines[0].startswith('@') and len(lines) > 3 and lines[2].startswith('+'):
        return 'fastq', ''
    if lines[0].startswith('>'):
        heads = [l for l in lines if l.startswith('>')]
        seqs = [l for l in lines if not l.startswith('>')]
        counted = sum(bool(re.search(r'(\d+)\s*$|[_x\-:|](\d+)$|count[=_:]?\d+', h)) for h in heads)
        if seqs and all(SEQ.match(s.strip()) for s in seqs[:50]):
            return ('fasta_collapsed' if counted >= 0.9 * len(heads) else 'fasta_reads'), heads[0][:60]
        return 'unknown', heads[0][:60]
    if re.match(r'^@(HD|SQ|PG)\b', lines[0]) or sum(len(l.split('\t')) >= 11 for l in lines) > 0.8 * len(lines):
        return 'sam_or_alignment', ''
    def is_data(line):
        f = [x for x in re.split(r'[\t,]', line.strip()) if x != '']
        return len(f) >= 2 and ((SEQ.match(f[0]) and re.match(r'^\d+(\.\d+)?$', f[1])) or
                                (re.match(r'^\d+(\.\d+)?$', f[0]) and SEQ.match(f[1])))
    body = lines if is_data(lines[0]) or not is_data(lines[1] if len(lines) > 1 else '') else lines[1:]
    sc = cs = named = 0
    for l in body[:150]:
        f = [x for x in re.split(r'[\t,]', l.strip()) if x != '']
        if len(f) >= 2 and SEQ.match(f[0]) and re.match(r'^\d+(\.\d+)?$', f[1]):
            sc += 1
        elif len(f) >= 2 and re.match(r'^\d+(\.\d+)?$', f[0]) and SEQ.match(f[1]):
            cs += 1
        elif f and not SEQ.match(f[0]):
            named += 1
    n = max(1, min(len(body), 150))
    header = '' if body is lines else lines[0][:80]
    if sc >= 0.9 * n:
        return 'tagcount_seq_count', header
    if cs >= 0.9 * n:
        return 'tagcount_count_seq', header
    if named >= 0.9 * n:
        return 'table_named', header
    return 'unknown', header


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', required=True, help='annotated sample TSV')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.samples), delimiter='\t'))
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\')
        w.writerow(['gse', 'gsm', 'file', 'format', 'header_or_first_id', 'n_columns', 'max_len', 'min_len',
                    'has_U', 'first_lines'])
        for r in rows:
            for url in [f for f in r['supplementary_files'].split(';') if f and f != 'NONE']:
                name = url.split('/')[-1]
                if not re.search(r'\.(txt|tsv|csv|fa|fasta|fna|fq|fastq|tag|seq|counts?)(\.gz)?$', name, re.I):
                    continue
                try:
                    text = decompress_partial(head_bytes(url), name.lower()).decode('utf-8', 'replace')
                except Exception as e:  # noqa: BLE001
                    w.writerow([r['gse'], r['gsm'], name, f'error:{type(e).__name__}', '', '', '', '', '', ''])
                    continue
                lines = text.split('\n')[:-1]
                fmt, hdr = classify(lines)
                seqs = [t for l in lines[:2000] for t in re.split(r'[\t, ]+', l.strip()) if SEQ.match(t)]
                ncol = len(re.split(r'[\t,]', lines[1].strip())) if len(lines) > 1 else ''
                w.writerow([r['gse'], r['gsm'], name, fmt, hdr, ncol, max(map(len, seqs)) if seqs else '',
                            min(map(len, seqs)) if seqs else '', any('U' in s.upper() for s in seqs[:500]),
                            ' || '.join(l[:70].replace('\t', '\\t').replace('"', "'") for l in lines[:4])])
                fh.flush()
                print(r['gsm'], name, fmt, file=sys.stderr)


if __name__ == '__main__':
    main()
