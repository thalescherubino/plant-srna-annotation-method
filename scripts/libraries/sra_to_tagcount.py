#!/usr/bin/env python3
"""Raw SRA small-RNA experiment -> trimmed, length-filtered Phasis tag-count file, with evidence-based trimming.

Added 2026-11-09 for tomato (user decision: build all libraries from raw SRA with one uniform pipeline, and take
care with kit-specific artefacts: 3' adapters, extra 5' bases, random 4N ends, inline barcodes, UMIs, poly-A
tailing; expect 21- and 24-nt peaks).

Per experiment (SRX, one library; all its runs; read 1 only for paired runs):
 1. download: prefetch + fasterq-dump, read 1 kept as <raw>/<SRX>.R1.fastq.gz (provenance, re-trimming)
 2. detect: on a 200,000-read sample, measure known 3' adapters (k-mer counts), a de novo over-represented
    12-mer when none is known, poly-A tails, and an invariant 5' prefix (inline barcode / ligation base).
 3. choose: evaluate every combination of {adapter or none} x {5' cut 0-4, prefix length} x {extra 3' cut 0-4} with
    cutadapt on the sample, counting reads of 18-34 nt that map to the genome with zero mismatches (Bowtie 1
    -v 0 -k 1). Pick the adapter with the most mapped reads; then, starting from no extra cut, accept a further
    cut only when it raises the mapped count by >= 20% (a real artefact: leftover 4N bases, barcodes, a fixed
    ligation base make reads unmappable, so removing them gives a large jump), taking the smallest cut within 3%
    of the best option; repeat until no cut gives such a jump. Small gains are NOT taken: cutting one more base
    always rescues a few reads with 3' untemplated additions or sequencing errors, and following those gains
    shifts the true 21/22/24-nt peaks down (tested 2026-11-09 on SRX3793170 and ERX1066761).
 4. trim the full library with the chosen settings (reads without the adapter are discarded when an adapter is
    used: they are inserts longer than the read), keep 18-34 nt, no N, collapse to SEQ<TAB>COUNT (DNA, no header).
 5. report: chosen settings, the whole evaluation table, length histogram, 5' nucleotide by length, mapping rate
    of the chosen sample, and FLAGS for review (no 21/24 peak, low mapping, few reads kept, UMI design, etc.).
Writes <out>/<SRX>.tag, <out>/<SRX>.report.json and one TSV row of statistics (stdout or --stats-row).

Usage: sra_to_tagcount.py --srx SRX... --runs SRR1;SRR2 --raw DIR --out DIR --index BOWTIE_INDEX [--threads 6]
Needs sra-tools, cutadapt, bowtie (1.x) and pigz on PATH (env plant-srna-reads).
"""
import argparse
import collections
import gzip
import hashlib
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

KNOWN_ADAPTERS = {
    'truseq_small_rna_RA3': 'TGGAATTCTCGGGTGCCAAGG',
    'illumina_small_rna_v1.5': 'ATCTCGTATGCCGTCTTCTGCTTG',
    'truseq_universal_nebnext': 'AGATCGGAAGAGCACACGTCT',
    'qiaseq_mirna': 'AACTGTAGGCACCATCAAT',
    'bgi_mgi_small_rna': 'AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA',
    'solid_ion_like': 'CGCCTTGGCCGTACAGCAG',
}
SAMPLE_N = 100_000
MINLEN, MAXLEN = 18, 34


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def download(srx, runs, raw_dir, threads, tmp):
    out = os.path.join(raw_dir, f'{srx}.R1.fastq.gz')
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out, 'existing'
    parts = []
    for srr in runs:
        d = os.path.join(tmp, srr)
        os.makedirs(d, exist_ok=True)
        run(['prefetch', '--max-size', '200G', '-O', d, srr], stdout=subprocess.DEVNULL)
        sra = os.path.join(d, srr, srr + '.sra')
        src = sra if os.path.exists(sra) else srr
        run(['fasterq-dump', '--split-files', '--skip-technical', '-e', str(threads), '-t', d, '-O', d, src],
            stdout=subprocess.DEVNULL)
        r1 = [os.path.join(d, f) for f in (f'{srr}_1.fastq', f'{srr}.fastq') if os.path.exists(os.path.join(d, f))]
        if not r1:
            raise RuntimeError(f'no read-1 fastq for {srr}')
        parts.append(r1[0])
    tmp_out = out + '.tmp'
    with open(tmp_out, 'wb') as fh:
        cat = subprocess.Popen(['cat'] + parts, stdout=subprocess.PIPE)
        run(['pigz', '-p', str(threads)], stdin=cat.stdout, stdout=fh)
        cat.wait()
    os.replace(tmp_out, out)
    for srr in runs:
        shutil.rmtree(os.path.join(tmp, srr), ignore_errors=True)
    return out, 'downloaded'


def read_sample(fq, n, skip=20_000):
    seqs = []
    with gzip.open(fq, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                k = i // 4
                if k < skip:
                    continue
                seqs.append(line.strip())
                if len(seqs) >= n:
                    break
    if len(seqs) < n // 2:  # small library: take from the start
        seqs = []
        with gzip.open(fq, 'rt') as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    seqs.append(line.strip())
                    if len(seqs) >= n:
                        break
    return seqs


def detect(seqs):
    det = {'n_sample': len(seqs), 'read_len_mode': collections.Counter(map(len, seqs)).most_common(1)[0][0]}
    known = {}
    for name, ad in KNOWN_ADAPTERS.items():
        k = ad[:10]
        known[name] = round(sum(1 for s in seqs if k in s) / len(seqs), 4)
    det['known_adapter_fraction'] = known
    det['polyA_fraction'] = round(sum(1 for s in seqs if 'AAAAAAAAAA' in s[15:]) / len(seqs), 4)
    # de novo: most common 12-mer starting at positions 16-40 (where a 3' adapter begins after an sRNA insert)
    c = collections.Counter()
    for s in seqs[:50_000]:
        for p in range(16, min(len(s) - 12, 40) + 1):
            c[s[p:p + 12]] += 1
    c = collections.Counter({k: v for k, v in c.items() if 'N' not in k and len(set(k)) > 2})  # no N or low-complexity runs
    top, cnt = c.most_common(1)[0] if c else ('', 0)
    det['denovo_12mer'], det['denovo_12mer_reads_frac'] = top, round(cnt / max(1, min(len(seqs), 50_000)), 4)
    # invariant 5' prefix: consecutive leading positions where one base > 85%
    L = 0
    for p in range(12):
        col = collections.Counter(s[p] for s in seqs if len(s) > p)
        b, n = col.most_common(1)[0]
        if n / sum(col.values()) > 0.85:
            L += 1
        else:
            break
    det['invariant_5prime_prefix_len'] = L
    det['first_base_composition'] = {b: round(v / len(seqs), 3) for b, v in collections.Counter(s[0] for s in seqs).items()}
    return det


def candidate_settings(det):
    ads = [(n, KNOWN_ADAPTERS[n]) for n, f in det['known_adapter_fraction'].items() if f >= 0.05]
    if not ads and det['denovo_12mer_reads_frac'] >= 0.05:
        ads = [('denovo_' + det['denovo_12mer'], det['denovo_12mer'])]
    if det['polyA_fraction'] >= 0.2:
        ads.append(('polyA', 'A' * 12))
    ads.append(('none', ''))
    cut5 = sorted({0, 1, 2, 3, 4, det['invariant_5prime_prefix_len']})
    cut5 = [c for c in cut5 if c <= 12]
    return [(a, c5, c3) for a, c5, c3 in itertools.product(ads, cut5, (0, 1, 2, 3, 4))]


def trim_two_pass(adapter, cut5, cut3, threads, inp, out, tmp):
    """adapter removal (+5' cut) first, then the extra 3' cut, then the length filter."""
    mid = os.path.join(tmp, 'mid.fastq.gz' if out.endswith('.gz') else 'mid.fastq')
    c1 = ['cutadapt', '-j', str(threads), '--quiet', '--max-n', '0', '-m', '1']
    if adapter:
        c1 += ['-a', adapter, '-e', '0.1', '-O', '6', '--discard-untrimmed']
    if cut5:
        c1 += ['-u', str(cut5)]
    run(c1 + ['-o', mid, inp])
    c2 = ['cutadapt', '-j', str(threads), '--quiet', '-m', str(MINLEN), '-M', str(MAXLEN)]
    if cut3:
        c2 += ['-u', str(-cut3)]
    run(c2 + ['-o', out, mid])
    os.remove(mid)


def evaluate(seqs, settings, index, threads, tmp):
    fq = os.path.join(tmp, 'sample.fastq')
    with open(fq, 'w') as fh:
        for i, s in enumerate(seqs):
            fh.write(f'@r{i}\n{s}\n+\n{"I" * len(s)}\n')
    rows = []
    for (aname, ad), c5, c3 in settings:
        out = os.path.join(tmp, 'eval.fastq')
        trim_two_pass(ad, c5, c3, 1, fq, out, tmp)
        lens, reads = collections.Counter(), []
        with open(out) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    s = line.strip()
                    lens[len(s)] += 1
                    reads.append(s)
        kept = len(reads)
        mapped = collections.Counter()
        if kept:
            fa = os.path.join(tmp, 'eval.fa')
            with open(fa, 'w') as fh:
                for i, s in enumerate(reads):
                    fh.write(f'>{i}\n{s}\n')
            p = subprocess.run(['bowtie', '-f', '-v', '0', '-k', '1', '-p', str(threads), '--suppress', '2,3,4,5,6,7,8',
                                index, fa], capture_output=True, text=True, check=True)
            for line in p.stdout.splitlines():
                mapped[len(reads[int(line.split('\t')[0])])] += 1
        m2024 = sum(mapped[L] for L in range(20, 25))
        rows.append({'adapter': aname, 'adapter_seq': ad, 'cut5': c5, 'cut3': c3, 'kept': kept,
                     'kept_frac': round(kept / len(seqs), 4), 'mapped_total': sum(mapped.values()),
                     'map_rate_total': round(sum(mapped.values()) / max(1, kept), 4), 'mapped_20_24': m2024,
                     'map_rate_20_24': round(m2024 / max(1, sum(lens[L] for L in range(20, 25))), 4),
                     'len21': lens[21], 'len24': lens[24],
                     'peak_len': lens.most_common(1)[0][0] if lens else 0})
    return rows


JUMP, NEAR = 1.20, 0.97


def choose(rows):
    """Adapter with the most mapped reads; then stepwise cuts, each accepted only for a >= 20% jump."""
    best_by_ad = collections.defaultdict(int)
    for r in rows:
        best_by_ad[r['adapter']] = max(best_by_ad[r['adapter']], r['mapped_total'])
    ad = max(best_by_ad, key=best_by_ad.get)
    grid = {(r['cut5'], r['cut3']): r for r in rows if r['adapter'] == ad}
    cur, path = (0, 0), []
    while True:
        here = grid[cur]['mapped_total']
        opts = [r for (c5, c3), r in grid.items() if c5 >= cur[0] and c3 >= cur[1] and (c5, c3) != cur]
        if not opts:
            break
        top = max(r['mapped_total'] for r in opts)
        if top < JUMP * max(1, here):
            break
        ok = [r for r in opts if r['mapped_total'] >= NEAR * top]
        nxt = min(ok, key=lambda r: (r['cut5'] + r['cut3'] - cur[0] - cur[1], -r['mapped_total']))
        path.append(f"{cur}->{(nxt['cut5'], nxt['cut3'])}: {here}->{nxt['mapped_total']}")
        cur = (nxt['cut5'], nxt['cut3'])
    return grid[cur], path


def collapse_to_tag(trimmed_fq, tag_out):
    c = collections.Counter()
    with gzip.open(trimmed_fq, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                c[line.strip().upper()] += 1
    with open(tag_out + '.tmp', 'w') as o:
        for s, n in sorted(c.items(), key=lambda x: (-x[1], x[0])):
            if re.fullmatch(r'[ACGT]+', s):
                o.write(f'{s}\t{n}\n')
    os.replace(tag_out + '.tmp', tag_out)
    return c


def stats_from_tag(tag):
    lens, five = collections.Counter(), collections.defaultdict(collections.Counter)
    total = distinct = 0
    with open(tag) as fh:
        for line in fh:
            s, n = line.rstrip('\n').split('\t')
            n = int(n)
            total += n
            distinct += 1
            lens[len(s)] += n
            five[len(s)][s[0]] += n
    f = lambda L: round(lens[L] / total, 4) if total else 0  # noqa: E731
    fp = lambda L: ','.join(f'{b}:{five[L][b] / max(1, sum(five[L].values())):.2f}' for b in 'ACGT')  # noqa: E731
    return total, distinct, lens, f, fp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--srx', required=True)
    ap.add_argument('--runs', required=True, help="';'-separated run accessions")
    ap.add_argument('--raw', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--index', required=True, help='Bowtie 1 genome index prefix')
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--tmp', default=None)
    a = ap.parse_args()
    os.makedirs(a.raw, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix=a.srx + '_', dir=a.tmp or a.out)
    rep = {'srx': a.srx, 'runs': a.runs}
    try:
        fq, how = download(a.srx, [r for r in a.runs.split(';') if r], a.raw, a.threads, tmp)
        rep['raw_fastq'], rep['download'] = fq, how
        seqs = read_sample(fq, SAMPLE_N)
        rep['detect'] = detect(seqs)
        evals = evaluate(seqs, candidate_settings(rep['detect']), a.index, a.threads, tmp)
        rep['evaluation'] = evals
        best, rep['choice_path'] = choose(evals)
        rep['chosen'] = {k: best[k] for k in ('adapter', 'adapter_seq', 'cut5', 'cut3')}
        trimmed = os.path.join(tmp, 'trimmed.fastq.gz')
        trim_two_pass(best['adapter_seq'], best['cut5'], best['cut3'], a.threads, fq, trimmed, tmp)
        tag = os.path.join(a.out, f'{a.srx}.tag')
        collapse_to_tag(trimmed, tag)
        total, distinct, lens, f, fp = stats_from_tag(tag)
        raw_reads = int(subprocess.run(f'pigz -dc -p 2 {fq} | wc -l', shell=True, capture_output=True, text=True,
                                       check=True).stdout) // 4
        flags = []
        hist = {L: lens[L] for L in range(MINLEN, MAXLEN + 1)}
        peak = max(hist, key=hist.get) if total else 0
        local = lambda L: hist[L] >= max(hist[L - 1], hist[L + 1]) if MINLEN < L < MAXLEN else False  # noqa: E731
        if not (local(21) or local(22)):
            flags.append('no_21_22_peak')
        if not (local(24) or local(23)):
            flags.append('no_23_24_peak')
        if peak not in (21, 22, 24):
            flags.append(f'main_peak_{peak}nt')
        if best['map_rate_20_24'] < 0.5:
            flags.append(f"low_exact_map_20_24:{best['map_rate_20_24']}")
        if total < 1_000_000:
            flags.append(f'few_reads_kept:{total}')
        if raw_reads and total / raw_reads < 0.2:
            flags.append(f'low_keep_fraction:{total / raw_reads:.2f}')
        rivals = [r for r in evals if r['adapter'] != best['adapter'] and r['mapped_total'] >= 0.9 * best['mapped_total']]
        if rivals:
            flags.append('ambiguous_adapter:' + rivals[0]['adapter'])
        if best['adapter'] == 'none' and rep['detect']['read_len_mode'] > 36:
            flags.append('no_adapter_found_long_reads')
        rep['flags'] = flags
        sha = hashlib.sha256(open(tag, 'rb').read()).hexdigest()
        row = {'gsm': a.srx, 'layout': f"sra_fastq_R1;adapter={best['adapter']};cut5={best['cut5']};cut3={best['cut3']}",
               'integer_counts': 'True', 'tagcount_file': tag, 'tagcount_sha256': sha, 'raw_reads': raw_reads,
               'kept_reads': total, 'kept_distinct': distinct, 'min_len': min(lens) if lens else 0,
               'max_len': max(lens) if lens else 0, 'frac_18_30': round(sum(lens[L] for L in range(18, 31)) / max(1, total), 4),
               'frac_21': f(21), 'frac_22': f(22), 'frac_24': f(24), 'five_prime_21': fp(21), 'five_prime_24': fp(24),
               'map_rate_20_24_sample': best['map_rate_20_24'], 'map_rate_total_sample': best['map_rate_total'],
               'flags': ';'.join(flags)}
        rep['stats'] = row
        rep['length_hist'] = hist
        json.dump(rep, open(os.path.join(a.out, f'{a.srx}.report.json'), 'w'), indent=1)
        print('\t'.join(str(row[k]) for k in row))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()
