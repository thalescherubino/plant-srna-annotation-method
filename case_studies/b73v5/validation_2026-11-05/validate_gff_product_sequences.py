#!/usr/bin/env python3
"""Validate GFF child feature coordinates against a local reference and class FASTA.

Requirements: Python 3 and bedtools. Genome must be an uncompressed indexed FASTA
with chromosome IDs matching the GFF. Does not edit any input files.

CAUTION: This is a class-level sequence-at-coordinate check, not an ID-aware
feature-to-FASTA association check. Match in either orientation is expected if a
reference-forward sRNA export was used; it does not by itself validate biology.
"""
import argparse
import collections
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path

RC = str.maketrans('ACGTN', 'TGCAN')

def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')

def read_product_sequences(path):
    sequences = set()
    parts = []
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if parts:
                    sequences.add(''.join(parts).upper())
                    count += 1
                    parts = []
            else:
                parts.append(''.join(line.split()))
        if parts:
            sequences.add(''.join(parts).upper())
            count += 1
    return sequences, count

def feature_id(attributes, default):
    for field in attributes.split(';'):
        if field.startswith('ID='):
            return field[3:]
    return default

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', required=True, help='Uncompressed FASTA indexed for bedtools; do not change headers')
    parser.add_argument('--gff', required=True, help='Curated GFF3 (1-based, closed)')
    parser.add_argument('--products', required=True, help='Class-specific canonical product FASTA (.fa or .fa.gz)')
    parser.add_argument('--feature-type', default='exon', help='GFF3 feature type to check (default: exon)')
    parser.add_argument('--out-prefix', required=True)
    args = parser.parse_args()
    if args.genome.endswith('.gz'):
        parser.error('Use an uncompressed genome FASTA for bedtools getfasta.')
    if not shutil_which('bedtools'):
        parser.error('bedtools is not available on PATH')

    products, n_products = read_product_sequences(args.products)
    details = {}
    counts = collections.Counter()
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    detail_path = str(prefix) + '.details.tsv'
    summary_path = str(prefix) + '.summary.txt'

    with tempfile.TemporaryDirectory(prefix='gff_coords_') as tempdir:
        bedfile = os.path.join(tempdir, 'features.bed')
        extracted = os.path.join(tempdir, 'extracted.bedseq')
        with open_text(args.gff) as source, open(bedfile, 'w') as out:
            for lineno, line in enumerate(source, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 9:
                    raise ValueError(f'Malformed GFF3 line {lineno}: expected 9 fields')
                if fields[2] != args.feature_type:
                    continue
                chrom, _, _, start, end, _, strand, _, attrs = fields
                start, end = int(start), int(end)
                if start < 1 or end < start:
                    raise ValueError(f'Bad GFF coordinates at line {lineno}: {start}-{end}')
                key = f'F{len(details)+1:08d}'
                details[key] = (feature_id(attrs, key), chrom, start, end, strand)
                print(chrom, start - 1, end, key, '0', strand, sep='\t', file=out)

        if not details:
            raise ValueError('No GFF features selected. Inspect the GFF feature types.')

        with open(extracted, 'w') as out:
            subprocess.run(['bedtools', 'getfasta', '-fi', args.genome,
                            '-bed', bedfile, '-bedOut'], stdout=out, check=True)

        seen = set()
        with open(detail_path, 'w') as out:
            header = ('feature_id', 'chrom', 'start_1based', 'end_1based', 'gff_strand',
                      'length', 'status', 'genome_forward_seq', 'genome_revcomp_seq')
            print('\t'.join(header), file=out)
            with open(extracted) as source:
                for lineno, line in enumerate(source, 1):
                    f = line.rstrip('\n').split('\t')
                    if len(f) < 7:
                        raise ValueError(f'Unexpected bedtools -bedOut row {lineno}')
                    key = f[3]
                    if key not in details or key in seen:
                        raise ValueError(f'Invalid/duplicate extracted feature key {key}')
                    seen.add(key)
                    fid, chrom, start, end, strand = details[key]
                    seq = f[-1].upper()
                    rc = seq.translate(RC)[::-1]
                    if len(seq) != end - start + 1:
                        status = 'EXTRACTION_LENGTH_MISMATCH'
                    else:
                        fwd, rev = seq in products, rc in products
                        status = ('BOTH_ORIENTATIONS' if fwd and rev else
                                  'FORWARD_ONLY' if fwd else
                                  'REVERSE_COMPLEMENT_ONLY' if rev else
                                  'NOT_FOUND_IN_PRODUCT_FASTA')
                    counts[(strand, status)] += 1
                    print(fid, chrom, start, end, strand, end - start + 1,
                          status, seq, rc, sep='\t', file=out)
            for key, (fid, chrom, start, end, strand) in details.items():
                if key not in seen:
                    counts[(strand, 'NO_REFERENCE_SEQUENCE')] += 1
                    print(fid, chrom, start, end, strand, end-start+1,
                          'NO_REFERENCE_SEQUENCE', '.', '.', sep='\t', file=out)

    ok = ('FORWARD_ONLY', 'REVERSE_COMPLEMENT_ONLY', 'BOTH_ORIENTATIONS')
    with open(summary_path, 'w') as out:
        print('GFF feature type:', args.feature_type, file=out)
        print('Selected GFF features:', len(details), file=out)
        print('Product FASTA records:', n_products, file=out)
        print('Unique product sequences:', len(products), file=out)
        matched = sum(v for (strand, status), v in counts.items() if status in ok)
        print('GFF intervals with a forward/RC sequence found in product FASTA:', matched, file=out)
        print('Percent compatible: %.4f%%' % (100 * matched / len(details)), file=out)
        print('IMPORTANT: Class-level sequence membership only. Does not prove each', file=out)
        print('GFF feature is associated with its exact FASTA query ID.', file=out)
        print('strand\tstatus\tcount', file=out)
        for (strand, status), value in sorted(counts.items()):
            print(f'{strand}\t{status}\t{value}', file=out)
    with open(summary_path) as result:
        print(result.read(), end='')
    print('Details:', detail_path)
    print('Summary:', summary_path)


def shutil_which(name):
    import shutil
    return shutil.which(name)

if __name__ == '__main__':
    main()
