#!/usr/bin/env python3
"""Inventory raw small-RNA sequencing experiments for one species in NCBI SRA (added 2026-11-09, tomato).

Used when processed GEO files are insufficient (tomato: 21 read-level GEO samples; user decision 2026-11-09:
build every library from raw SRA reads with one uniform pipeline). One output row per SRA experiment (SRX =
one library; its runs are concatenated later). Columns follow geo_inventory.py, so curate_geo_candidates.py and
select_libraries.py work unchanged:
  gsm      = SRX experiment accession (the library id)      gse = SRP study accession (per-series cap)
  title, source, characteristics (BioSample attributes "tag: value | ..."), organism, library_strategy,
  library_selection, supplementary_files (= the run accessions, ';'-joined)
plus SRA-specific columns: runs, layout, platform, instrument, spots, bases, avg_len, design, protocol,
geo_sample (GSM alias if the experiment came from GEO), sample_accession, bioproject.
The design/protocol text is kept because it often names the library kit (adapters, UMIs, random bases).

Usage: sra_inventory.py --organism "Solanum lycopersicum" --out inventory/sra
       [--term '<custom SRA query>']   (default: organism AND (miRNA-Seq OR ncRNA-Seq strategy))
Writes <out>/sra_query.json, <out>/sra_xml/*.xml (raw records) and <out>/geo_samples.tsv.
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def get(url, data=None, tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'User-Agent': 'plant-srna-annotation-method/sra_inventory'})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                raise
            print(f'retry {k + 1} after {type(e).__name__}', file=sys.stderr)
            time.sleep(5 * (k + 1))


def eutil(name, **params):
    time.sleep(0.4)  # stay under the 3 requests/s anonymous E-utilities limit
    return get(EUTILS + name + '.fcgi', data=urllib.parse.urlencode(params).encode())


def text(el, path):
    x = el.find(path)
    return (x.text or '').strip() if x is not None and x.text else ''


def parse_package(p):
    exp, study, sample = p.find('EXPERIMENT'), p.find('STUDY'), p.find('SAMPLE')
    lib = exp.find('DESIGN/LIBRARY_DESCRIPTOR')
    layout = lib.find('LIBRARY_LAYOUT')
    plat = exp.find('PLATFORM')
    platform = plat[0].tag if plat is not None and len(plat) else ''
    attrs = []
    for a in sample.findall('SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE') if sample is not None else []:
        attrs.append(f"{text(a, 'TAG')}: {text(a, 'VALUE')}")
    geo = ''
    for x in (sample.findall('IDENTIFIERS/EXTERNAL_ID') if sample is not None else []):
        if (x.get('namespace') or '').upper() == 'GEO':
            geo = (x.text or '').strip()
    if not geo and sample is not None and (sample.get('alias') or '').startswith('GSM'):
        geo = sample.get('alias')
    runs = p.findall('RUN_SET/RUN')
    spots = sum(int(r.get('total_spots') or 0) for r in runs)
    bases = sum(int(r.get('total_bases') or 0) for r in runs)
    proj = ''
    for x in (study.findall('IDENTIFIERS/EXTERNAL_ID') if study is not None else []):
        if (x.get('namespace') or '').lower() == 'bioproject':
            proj = (x.text or '').strip()
    return {
        'gse': study.get('accession') if study is not None else '',
        'gsm': exp.get('accession'),
        'title': exp.findtext('TITLE', '').strip(),
        'source': (text(sample, 'TITLE') if sample is not None else ''),
        'organism': text(sample, 'SAMPLE_NAME/SCIENTIFIC_NAME') if sample is not None else '',
        'characteristics': ' | '.join(attrs),
        'molecule': lib.findtext('LIBRARY_SOURCE', ''),
        'library_strategy': lib.findtext('LIBRARY_STRATEGY', ''),
        'library_selection': lib.findtext('LIBRARY_SELECTION', ''),
        'instrument': text(plat[0], 'INSTRUMENT_MODEL') if platform else '',
        'treatment': '',
        'supplementary_files': ';'.join(r.get('accession') for r in runs),
        'relation': '',
        'runs': ';'.join(r.get('accession') for r in runs),
        'layout': layout[0].tag if layout is not None and len(layout) else '',
        'platform': platform,
        'spots': spots,
        'bases': bases,
        'avg_len': round(bases / spots, 1) if spots else '',
        'design': text(exp, 'DESIGN/DESIGN_DESCRIPTION'),
        'protocol': lib.findtext('LIBRARY_CONSTRUCTION_PROTOCOL', '').strip(),
        'geo_sample': geo,
        'sample_accession': sample.get('accession') if sample is not None else '',
        'bioproject': proj,
        'study_title': text(study, 'DESCRIPTOR/STUDY_TITLE') if study is not None else '',
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--organism', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--term')
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, 'sra_xml'), exist_ok=True)
    term = a.term or f'"{a.organism}"[Organism] AND ("mirna seq"[Strategy] OR "ncrna seq"[Strategy])'
    res = json.loads(eutil('esearch', db='sra', term=term, retmax=10000, retmode='json'))
    ids = res['esearchresult']['idlist']
    json.dump({'term': term, 'count': res['esearchresult']['count'], 'ids': ids,
               'date': time.strftime('%Y-%m-%dT%H:%M:%S%z')}, open(os.path.join(a.out, 'sra_query.json'), 'w'), indent=1)
    rows = []
    for k in range(0, len(ids), 200):
        f = os.path.join(a.out, 'sra_xml', f'batch_{k // 200:03d}.xml')
        if not os.path.exists(f):
            open(f, 'wb').write(eutil('efetch', db='sra', id=','.join(ids[k:k + 200]), retmode='xml'))
        for p in ET.parse(f).getroot().findall('EXPERIMENT_PACKAGE'):
            rows.append(parse_package(p))
        print(f'batch {k // 200}: {len(rows)} experiments', file=sys.stderr)
    with open(os.path.join(a.out, 'geo_samples.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'experiments={len(rows)} studies={len({r["gse"] for r in rows})}', file=sys.stderr)


if __name__ == '__main__':
    main()
