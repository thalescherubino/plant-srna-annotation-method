#!/usr/bin/env python3
"""Inventory NCBI GEO small-RNA series for a species (read-only web queries).

1. esearch (db=gds) for series of type "non coding rna profiling by high throughput
   sequencing" for the organism;
2. esummary for series metadata (samples, platforms, supplementary file types);
3. series_matrix download per series -> one row per sample with title, source,
   characteristics, organism, molecule and supplementary files.

Respects NCBI's 3 requests/s limit without an API key. Writes:
  <out>/geo_series.tsv, <out>/geo_samples.tsv, <out>/geo_query.json, raw matrices.
"""
import argparse
import csv
import gzip
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
UA = {'User-Agent': 'plant-srna-annotation-inventory/0.1'}


def get(url, tries=4):
    for k in range(tries):
        try:
            time.sleep(0.4)
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


def eutil(name, **params):
    params.setdefault('tool', 'plant_srna_inventory')
    return json.loads(get(EUTILS + name + '.fcgi?' + urllib.parse.urlencode(params)))


def matrix_urls(gse):
    stub = gse[:-3] + 'nnn'
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/'
    html = get(base).decode('utf-8', 'replace')
    names = sorted({p.split('"')[0] for p in html.split('href="')[1:] if p.split('"')[0].endswith('_series_matrix.txt.gz')})
    return [base + n for n in names]


def parse_matrix(raw):
    rows = {}
    for line in io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(raw)), encoding='utf-8', errors='replace'):
        if not line.startswith('!Sample_'):
            continue
        parts = next(csv.reader([line.rstrip('\n')], delimiter='\t'))
        rows.setdefault(parts[0], []).append(parts[1:])
    gsms = rows.get('!Sample_geo_accession', [[]])[0]
    out = []
    for i, gsm in enumerate(gsms):
        def field(key, join=' | '):
            vals = [r[i] for r in rows.get(key, []) if i < len(r) and r[i]]
            return join.join(vals)
        out.append({
            'gsm': gsm, 'title': field('!Sample_title'), 'source': field('!Sample_source_name_ch1'),
            'organism': field('!Sample_organism_ch1'), 'characteristics': field('!Sample_characteristics_ch1'),
            'molecule': field('!Sample_molecule_ch1'), 'library_strategy': field('!Sample_library_strategy'),
            'library_selection': field('!Sample_library_selection'), 'instrument': field('!Sample_instrument_model'),
            'treatment': field('!Sample_treatment_protocol_ch1')[:300],
            'supplementary_files': ';'.join(v for k in sorted(rows) if k.startswith('!Sample_supplementary_file')
                                            for v in [field(k, ';')] if v),
            'relation': field('!Sample_relation', ';'),
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--organism', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--term', help='custom GEO DataSets query (default: organism + ncRNA-seq dataset type)')
    ap.add_argument('--reparse', action='store_true',
                    help='rebuild geo_samples.tsv from already-downloaded matrices (no network)')
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, 'matrices'), exist_ok=True)
    if a.reparse:
        import glob
        series = {r['gse']: r for r in csv.DictReader(open(os.path.join(a.out, 'geo_series.tsv')), delimiter='\t')}
        sample_rows = []
        for p in sorted(glob.glob(os.path.join(a.out, 'matrices', '*_series_matrix.txt.gz'))):
            gse = os.path.basename(p).split('_')[0].split('-')[0]
            for r in parse_matrix(open(p, 'rb').read()):
                r['gse'] = gse
                sample_rows.append(r)
        write_samples(a.out, sample_rows)
        print(f'reparsed series={len(series)} samples={len(sample_rows)}', file=sys.stderr)
        return
    term = a.term or (f'"{a.organism}"[Organism] AND gse[ETYP] AND '
                      '"non coding rna profiling by high throughput sequencing"[DataSet Type]')
    res = eutil('esearch', db='gds', term=term, retmax=5000, retmode='json')
    ids = res['esearchresult']['idlist']
    json.dump({'term': term, 'count': res['esearchresult']['count'], 'ids': ids,
               'date': time.strftime('%Y-%m-%dT%H:%M:%S%z')}, open(os.path.join(a.out, 'geo_query.json'), 'w'), indent=1)
    series = []
    for k in range(0, len(ids), 100):
        summ = eutil('esummary', db='gds', id=','.join(ids[k:k + 100]), retmode='json')['result']
        for uid in summ.get('uids', []):
            s = summ[uid]
            series.append({'gse': 'GSE' + s['gse'], 'title': s['title'], 'n_samples': s['n_samples'],
                           'taxon': s['taxon'], 'gdstype': s['gdstype'], 'suppfile': s.get('suppfile', ''),
                           'pdat': s['pdat'], 'pubmed': ','.join(map(str, s.get('pubmedids', [])))})
    sample_rows = []
    for s in series:
        try:
            urls = matrix_urls(s['gse'])
        except Exception as e:  # noqa: BLE001
            s['matrix_status'] = f'error:{type(e).__name__}'
            continue
        n = 0
        for u in urls:
            raw = get(u)
            open(os.path.join(a.out, 'matrices', os.path.basename(u)), 'wb').write(raw)
            for r in parse_matrix(raw):
                r['gse'] = s['gse']
                sample_rows.append(r)
                n += 1
        s['matrix_status'] = f'ok:{len(urls)}_files:{n}_samples'
        print(s['gse'], s['n_samples'], s['matrix_status'], file=sys.stderr)
    with open(os.path.join(a.out, 'geo_series.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(series[0].keys()) + ['matrix_status'], delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(series)
    write_samples(a.out, sample_rows)
    print(f'series={len(series)} samples={len(sample_rows)}', file=sys.stderr)


def write_samples(out, sample_rows):
    fields = ['gse', 'gsm', 'title', 'source', 'organism', 'characteristics', 'molecule', 'library_strategy',
              'library_selection', 'instrument', 'treatment', 'supplementary_files', 'relation']
    with open(os.path.join(out, 'geo_samples.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(sample_rows)


if __name__ == '__main__':
    main()
