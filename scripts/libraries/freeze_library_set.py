#!/usr/bin/env python3
"""Freeze a selected library set: tag-count links, manifest with provenance/stats, ShortStack condensed FASTA
(headers >LIB_Cd<n>_<count>, the ShortStack 4.1 condensed format), and miRador library links.

Usage: freeze_library_set.py <selection.selected.tsv> <outdir> <libraries_dir> [--sra-inventory <annotated.tsv>]
GEO mode: <libraries_dir> must contain library_stats.tsv, tagcount/ and geo_raw/download_manifest.tsv
(outputs of inspect_and_standardize_libraries.py and download_geo_libraries.py).
SRA mode (--sra-inventory, added 2026-11-09): <libraries_dir> holds library_stats.tsv from sra_to_tagcount.py
(absolute tagcount paths); provenance is the SRA run accessions, the trimming chosen and the QC flags."""
import csv
import gzip
import hashlib
import os
import sys

sel_path, F, LIBDIR = sys.argv[1], sys.argv[2], sys.argv[3]
SRA = None
if '--sra-inventory' in sys.argv:
    SRA = {r['gsm']: r for r in csv.DictReader(open(sys.argv[sys.argv.index('--sra-inventory') + 1]), delimiter='\t')}
for sub in ('tag', 'shortstack_condensed', 'mirador_libs'):
    os.makedirs(os.path.join(F, sub), exist_ok=True)
sel = list(csv.DictReader(open(sel_path), delimiter='\t'))
stats = {s['gsm']: s for s in csv.DictReader(open(LIBDIR + '/library_stats.tsv'), delimiter='\t')}
man = {} if SRA else {m['gsm']: m for m in csv.DictReader(open(LIBDIR + '/geo_raw/download_manifest.tsv'), delimiter='\t')}
rows = []
for r in sel:
    s = stats[r['gsm']]
    src = s['tagcount_file'] if SRA else LIBDIR + '/tagcount/' + s['tagcount_file']
    lib = r['gsm']
    if hashlib.sha256(open(src, 'rb').read()).hexdigest() != s['tagcount_sha256']:
        sys.exit(f'checksum mismatch {src}')
    for link, target in ((f'{F}/tag/{lib}.tag', src), (f'{F}/mirador_libs/{lib}.txt', f'../tag/{lib}.tag')):
        if not os.path.lexists(link):
            os.symlink(target, link)
    cond = f'{F}/shortstack_condensed/{lib}.fa.gz'
    if not os.path.exists(cond):
        with open(src) as fin, gzip.open(cond + '.part', 'wt', compresslevel=1) as fout:
            for n, line in enumerate(fin, 1):
                seq, cnt = line.rstrip('\n').split('\t')
                fout.write(f'>{lib}_Cd{n}_{cnt}\n{seq}\n')
        os.replace(cond + '.part', cond)
    if SRA:
        a = SRA[lib]
        rows.append({'library_id': lib, 'selection_order': r['selection_order'], 'gse': r['gse'], 'title': r['title'],
                     'organ': r['organ'], 'subspecies': r['subspecies'], 'genotype_class': r['genotype_class'],
                     'treatment_class': r['treatment_class'], 'sra_runs': a['runs'], 'sra_layout': a['layout'],
                     'instrument': a['instrument'], 'study_title': a['study_title'], 'geo_sample': a['geo_sample'],
                     'trimming': s['layout'], 'raw_reads': s['raw_reads'], 'kept_reads': s['kept_reads'],
                     'kept_distinct': s['kept_distinct'], 'min_len': s['min_len'], 'max_len': s['max_len'],
                     'frac_21': s['frac_21'], 'frac_22': s['frac_22'], 'frac_24': s['frac_24'],
                     'five_prime_21': s['five_prime_21'], 'five_prime_24': s['five_prime_24'],
                     'map_rate_total_sample': s['map_rate_total_sample'], 'qc_flags': s['flags'],
                     'tagcount_file': f'tag/{lib}.tag', 'tagcount_sha256': s['tagcount_sha256'],
                     'characteristics': r['characteristics']})
        continue
    rows.append({'library_id': lib, 'selection_order': r['selection_order'], 'gse': r['gse'], 'title': r['title'],
                 'organ': r['organ'], 'subspecies': r['subspecies'], 'genotype_class': r['genotype_class'],
                 'treatment_class': r['treatment_class'], 'source_url': man[lib]['url'], 'source_sha256': man[lib]['sha256'],
                 'source_layout': s['layout'],
                 'conversion': 'none (already SEQ<TAB>COUNT)' if s['layout'] == 'seq_count' and not s['notes']
                 else f"standardized from {s['layout']}; {s['notes']}",
                 'prefiltered_flag': r['prefiltered_flag'], 'kept_reads': s['kept_reads'], 'kept_distinct': s['kept_distinct'],
                 'dropped_nonACGT_reads': s['dropped_nonACGT_reads'], 'dropped_length_reads': s['dropped_length_reads'],
                 'min_len': s['min_len'], 'max_len': s['max_len'], 'frac_21': s['frac_21'], 'frac_22': s['frac_22'],
                 'frac_24': s['frac_24'], 'five_prime_21': s['five_prime_21'], 'five_prime_24': s['five_prime_24'],
                 'tagcount_file': f'tag/{lib}.tag', 'tagcount_sha256': s['tagcount_sha256'],
                 'characteristics': r['characteristics']})
with open(F + '/library_manifest.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(len(rows), 'libraries;', sum(int(r['kept_reads']) for r in rows) / 1e6, 'M reads')
