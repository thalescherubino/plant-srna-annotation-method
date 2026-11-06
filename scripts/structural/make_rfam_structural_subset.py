#!/usr/bin/env python3
"""Build the Rfam structural-RNA model subset used for genome-wide structural RNA prediction.

Selects families whose Rfam type starts with "Gene; rRNA" or "Gene; snRNA" (spliceosomal snRNA, C/D-box and
H/ACA-box snoRNA, scaRNA); tRNA families are excluded because tRNAs are predicted with tRNAscan-SE.
For Rfam 15.1 this gives 827 families. Writes, in the Rfam release directory:
  structural_subset.families.tsv   rfam_id, accession, type
  structural_subset.acc.txt
  Rfam.structural_subset.cm (+ cmpress indexes)
  Rfam.structural_subset.clanin    Rfam.clanin restricted to subset members (clans with >= 2 members kept;
                                   cmscan --clanin requires every listed model to be present)
Requires Infernal (cmfetch, cmpress) on PATH and an uncompressed Rfam.cm with family.txt.gz in --rfam-dir.
"""
import argparse
import gzip
import os
import subprocess


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rfam-dir', required=True, help='directory with Rfam.cm, Rfam.clanin and family.txt.gz')
    a = ap.parse_args()
    d = a.rfam_dir
    fams = []
    with gzip.open(os.path.join(d, 'family.txt.gz'), 'rt', errors='replace') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) > 18 and (f[18].startswith('Gene; rRNA') or f[18].startswith('Gene; snRNA')):
                fams.append((f[1], f[0], f[18]))
    with open(os.path.join(d, 'structural_subset.families.tsv'), 'w') as out:
        out.writelines(f'{n}\t{acc}\t{t}\n' for n, acc, t in fams)
    with open(os.path.join(d, 'structural_subset.acc.txt'), 'w') as out:
        out.writelines(f'{acc}\n' for _, acc, _ in fams)
    cm = os.path.join(d, 'Rfam.cm')
    subprocess.run(['cmfetch', '--index', cm], check=True, stdout=subprocess.DEVNULL)
    sub = os.path.join(d, 'Rfam.structural_subset.cm')
    with open(sub, 'w') as out:
        subprocess.run(['cmfetch', '-f', cm, os.path.join(d, 'structural_subset.acc.txt')], check=True, stdout=out)
    subprocess.run(['cmpress', '-F', sub], check=True, stdout=subprocess.DEVNULL)
    keep = {n for n, _, _ in fams}
    n_clans = 0
    with open(os.path.join(d, 'Rfam.clanin')) as fin, open(os.path.join(d, 'Rfam.structural_subset.clanin'), 'w') as out:
        for line in fin:
            f = line.rstrip('\n').split('\t')
            members = [x for x in f[1:] if x in keep]
            if len(members) >= 2:
                out.write('\t'.join([f[0]] + members) + '\n')
                n_clans += 1
    print(f'families={len(fams)} clans_kept={n_clans}')


if __name__ == '__main__':
    main()
