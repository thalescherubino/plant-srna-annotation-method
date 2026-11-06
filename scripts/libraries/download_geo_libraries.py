#!/usr/bin/env python3
"""Download the selected GEO supplementary file of every eligible library (read-only source; resumable).

Usage: download_geo_libraries.py --pool <eligible_pool.tsv> --samples <geo_samples.annotated.tsv> --outdir <dir>
Writes <outdir>/<GSE>/<file> and appends <outdir>/download_manifest.tsv (gse, gsm, url, bytes, sha256, status).
"""
import argparse, csv, hashlib, os, sys, time, urllib.request

ap = argparse.ArgumentParser()
ap.add_argument('--pool', required=True, help='TSV with gse, gsm, selected_file columns')
ap.add_argument('--samples', required=True, help='annotated GEO sample table (supplementary_files column)')
ap.add_argument('--outdir', required=True)
a = ap.parse_args()
pool = list(csv.DictReader(open(a.pool), delimiter='\t'))
samples = {r['gsm']: r for r in csv.DictReader(open(a.samples), delimiter='\t')}
os.makedirs(a.outdir, exist_ok=True)
mpath = os.path.join(a.outdir, 'download_manifest.tsv')
if not os.path.exists(mpath):
    open(mpath, 'w').write('gse\tgsm\turl\tbytes\tsha256\tstatus\n')
done = {l.split('\t')[1] for l in open(mpath) if l.rstrip('\n').endswith('\tok')}
man = open(mpath, 'a')
for r in pool:
    if r['gsm'] in done:
        continue
    url = next(u for u in samples[r['gsm']]['supplementary_files'].split(';') if u.endswith('/' + r['selected_file']))
    url = url.replace('ftp://ftp.ncbi.nlm.nih.gov', 'https://ftp.ncbi.nlm.nih.gov')
    d = os.path.join(a.outdir, r['gse']); os.makedirs(d, exist_ok=True)
    path = os.path.join(d, r['selected_file'])
    status = 'ok'
    for k in range(5):
        try:
            time.sleep(0.5)
            req = urllib.request.Request(url, headers={'User-Agent': 'plant-srna-annotation-inventory/0.1'})
            with urllib.request.urlopen(req, timeout=600) as resp, open(path + '.part', 'wb') as fh:
                for b in iter(lambda: resp.read(1 << 20), b''):
                    fh.write(b)
            os.replace(path + '.part', path)
            break
        except Exception as e:  # noqa: BLE001
            status = f'error:{type(e).__name__}'
            time.sleep(5 * (k + 1))
    else:
        man.write(f"{r['gse']}\t{r['gsm']}\t{url}\t\t\t{status}\n"); man.flush(); continue
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    man.write(f"{r['gse']}\t{r['gsm']}\t{url}\t{os.path.getsize(path)}\t{h}\tok\n"); man.flush()
    print(r['gsm'], os.path.getsize(path), file=sys.stderr)
