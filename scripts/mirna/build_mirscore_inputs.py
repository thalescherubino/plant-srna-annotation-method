#!/usr/bin/env python3
"""Build miRScore hairpin/mature FASTA inputs for one miRNA candidate set.

miRScore pairs records by name: hairpin X, mature X.mature, star X.star (plant mode). Reads are mapped to the
hairpins only, so the hairpin span sets the precision denominator. For predictor sets the hairpin rule is the
same for every tool: duplex span (leftmost of miR/star to rightmost) +/- FLANK nt (default 20), cut from the
genome and reverse-complemented on the minus strand. FLANK is reduced so the hairpin never exceeds MAX_LEN
(default 300, the miRScore plant-mode precursor limit); a span already >= MAX_LEN gets no flank.

Sets:
  shortstack  ShortStack 4 results dir (Results.txt MIRNA=Y, mir.fasta .mature/.star coordinates)
  mirador     miRador finalAnnotatedCandidates.csv (Chr, Strand c=minus, 1-based leftmost miR/Star Position)
  mirbase     miRBase hairpin.fa + mature.fa, filtered to a species prefix, used as deposited. Hairpins whose
              mature/star touches the first/last base are extended from the genome when the hairpin is found
              exactly (miRScore exits on such records), using the same capped-flank rule; this is logged.
Writes <out>/hairpin.fa, <out>/mature.fa, <out>/id_map.tsv, <out>/build.log
"""
import argparse
import csv
import os
import re
import subprocess
import sys

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


class Fasta:
    def __init__(self, path):
        self.path = path
        self.idx = {}
        with open(path + '.fai') as fh:
            for line in fh:
                f = line.split('\t')
                self.idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        self.fh = open(path, 'rb')

    def fetch(self, chrom, start, end):  # 1-based closed
        ln, off, lb, lw = self.idx[chrom]
        start = max(1, start); end = min(ln, end)
        s0, e0 = start - 1, end
        b0 = off + (s0 // lb) * lw + s0 % lb
        b1 = off + (e0 // lb) * lw + e0 % lb
        self.fh.seek(b0)
        return self.fh.read(b1 - b0).replace(b'\n', b'').replace(b'\r', b'').decode().upper()

    def length(self, chrom):
        return self.idx[chrom][0]


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def dna(s):
    return s.upper().replace('U', 'T')


def capped_flanks(span, flank, max_len):
    """Flank each side by `flank` nt unless that exceeds max_len; then add only max_len - span nt in total
    (split evenly, the odd nt on the right). Returns (left, right); (0, 0) if span already >= max_len."""
    extra = max(0, min(2 * flank, max_len - span))
    return extra // 2, extra - extra // 2


def extend_locus(g, chrom, strand, lo, hi, flank, max_len):
    """lo..hi 1-based closed on + coordinates; returns transcript-oriented sequence and new lo, hi."""
    left, right = capped_flanks(hi - lo + 1, flank, max_len)
    if strand == '-':
        left, right = right, left  # keep the odd nt on the transcript 3' side
    lo = max(1, lo - left); hi = min(g.length(chrom), hi + right)
    s = g.fetch(chrom, lo, hi)
    return (rc(s) if strand == '-' else s), lo, hi


def duplex_hairpin(g, chrom, strand, a, b, flank, max_len):
    """a, b: (start, end) 1-based closed of miR and star on the + strand coordinates."""
    return extend_locus(g, chrom, strand, min(a[0], b[0]), max(a[1], b[1]), flank, max_len)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--set', required=True, choices=['shortstack', 'mirador', 'mirbase'])
    ap.add_argument('--genome', required=True)
    ap.add_argument('--input', required=True, help='ShortStack results dir | miRador CSV | miRBase dir')
    ap.add_argument('--species-prefix', default='osa-', help='miRBase species prefix')
    ap.add_argument('--flank', type=int, default=20, help='nt added on each side of the duplex (predictor sets)')
    ap.add_argument('--max-len', type=int, default=300, help='hairpin length cap; miRScore plant mode fails precursors > 300 nt')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    g = Fasta(a.genome)
    hp, mat, idmap, log = [], [], [], []

    def add(name, hseq, mseq, sseq, src, loc):
        if mseq not in hseq or (sseq and sseq not in hseq):
            log.append(f'SKIP {name} ({src}): mature/star not found in hairpin')
            return
        hp.append((name, hseq)); mat.append((name + '.mature', mseq))
        if sseq:
            mat.append((name + '.star', sseq))
        idmap.append((name, src, loc, len(hseq), mseq, sseq))

    if a.set == 'shortstack':
        keep = set()
        with open(os.path.join(a.input, 'Results.txt')) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['MIRNA'] == 'Y':
                    keep.add(r['Name'])
        coords = {}
        for h, _ in read_fasta(os.path.join(a.input, 'mir.fasta')):
            m = re.match(r'(\S+?)(\.mature|\.star)?::(\S+):(\d+)-(\d+)\(([+-])\)', h)
            if not m:
                sys.exit(f'unparsed header {h}')
            coords.setdefault(m.group(1), {})[m.group(2) or 'hp'] = (m.group(3), int(m.group(4)), int(m.group(5)), m.group(6))
        for name in sorted(keep, key=lambda x: int(x.split('_')[-1])):
            c = coords.get(name)
            if not c or '.mature' not in c or '.star' not in c:
                log.append(f'SKIP {name}: no mature/star in mir.fasta'); continue
            chrom, ms, me, strand = c['.mature']; _, ss, se, _ = c['.star']
            hseq, lo, hi = duplex_hairpin(g, chrom, strand, (ms, me), (ss, se), a.flank, a.max_len)
            mseq = g.fetch(chrom, ms, me); sseq = g.fetch(chrom, ss, se)
            if strand == '-':
                mseq, sseq = rc(mseq), rc(sseq)
            add(name, hseq, mseq, sseq, name, f'{chrom}:{lo}-{hi}({strand})')
        log.append(f'MIRNA=Y loci in Results.txt: {len(keep)}')

    elif a.set == 'mirador':
        with open(a.input) as fh:
            rows = list(csv.DictReader(fh))
        for i, r in enumerate(rows, 1):
            name = f'MRD{i:04d}'
            chrom = r['Chr']; strand = '-' if r['Strand'] in ('c', '-') else '+'
            mseq, sseq = dna(r['miR Sequence']), dna(r['Star Sequence'])
            ms, ss = int(r['miR Position']), int(r['Star Position'])
            mi = (ms, ms + len(mseq) - 1); si = (ss, ss + len(sseq) - 1)
            gm = g.fetch(chrom, *mi); gm = rc(gm) if strand == '-' else gm
            if gm != mseq:
                log.append(f'SKIP {name} ({r["miR Name"]}): genome at miR Position != miR Sequence'); continue
            hseq, lo, hi = duplex_hairpin(g, chrom, strand, mi, si, a.flank, a.max_len)
            add(name, hseq, mseq, sseq, r['miR Name'], f'{chrom}:{lo}-{hi}({strand})')
        log.append(f'miRador final candidates: {len(rows)}')

    else:  # mirbase
        pre = a.species_prefix
        hps = [(h.split()[0], dna(s)) for h, s in read_fasta(os.path.join(a.input, 'hairpin.fa')) if h.startswith(pre)]
        mats = [(h.split()[0], dna(s)) for h, s in read_fasta(os.path.join(a.input, 'mature.fa')) if h.startswith(pre)]
        # genome lookup for hairpins whose mature touches an end: exact search with bowtie is overkill for a
        # few records; scan chromosomes with str.find on both strands.
        chroms = {c: None for c in g.idx}
        mat_seqs = [s for _, s in mats]
        ext = 0
        for name, hseq in hps:
            touches = any(hseq.startswith(s) or hseq.endswith(s) for s in mat_seqs if s in hseq)
            if touches:
                loc = None
                for c in chroms:
                    if chroms[c] is None:
                        chroms[c] = g.fetch(c, 1, g.length(c))
                    for strand, q in (('+', hseq), ('-', rc(hseq))):
                        p = chroms[c].find(q)
                        if p >= 0:
                            loc = (c, p + 1, p + len(hseq), strand); break
                    if loc:
                        break
                if loc:
                    c, s, e, strand = loc
                    new, s2, e2 = extend_locus(g, c, strand, s, e, a.flank, a.max_len)
                    if len(new) == len(hseq):
                        log.append(f'NOT_EXTENDED {name}: mature/star at hairpin end but hairpin already {len(hseq)} nt (max {a.max_len})')
                    else:
                        log.append(f'EXTENDED {name}: mature/star at hairpin end; {len(hseq)} -> {len(new)} nt from {c}:{s}-{e}({strand})')
                        hseq = new; ext += 1
                else:
                    log.append(f'NOT_EXTENDED {name}: mature/star at hairpin end and hairpin not found exactly in genome')
            hp.append((name, hseq)); idmap.append((name, name, '', len(hseq), '', ''))
        mat = mats
        log.append(f'miRBase {pre} hairpins: {len(hps)}; matures: {len(mats)}; extended: {ext}')

    with open(os.path.join(a.out, 'hairpin.fa'), 'w') as out:
        out.writelines(f'>{n}\n{s}\n' for n, s in hp)
    with open(os.path.join(a.out, 'mature.fa'), 'w') as out:
        out.writelines(f'>{n}\n{s}\n' for n, s in mat)
    with open(os.path.join(a.out, 'id_map.tsv'), 'w') as out:
        out.write('mirscore_id\tsource_id\thairpin_locus\thairpin_len\tmature\tstar\n')
        out.writelines('\t'.join(map(str, r)) + '\n' for r in idmap)
    n_over = sum(len(s) > a.max_len for _, s in hp)
    log.append(f'written hairpins={len(hp)} mature_records={len(mat)} flank={a.flank} max_len={a.max_len} hairpins_over_max={n_over}')
    with open(os.path.join(a.out, 'build.log'), 'w') as out:
        out.writelines(l + '\n' for l in log)
    print('\n'.join(log[-3:]))


if __name__ == '__main__':
    main()
