#!/usr/bin/env python3
"""Annotate GEO samples of one species for small-RNA library selection (no downloads of data files).

Usage: curate_geo_candidates.py <inventory_dir> [--species rice|tomato]   (default rice, as originally written)
Species-specific rules (genotype groups, organ categories, extra mutant/special-line terms) are in species_profiles.py.
Adds derived columns to <inventory_dir>/geo_samples.tsv:
  subspecies      genotype group from the species profile (rice: japonica_nipponbare | japonica_other |
                  indica_or_aus | unknown_sativa; tomato: heinz_1706 | cultivar_other | unknown_lycopersicum |
                  cerasiforme | interspecific_line), or non_target
  srna_library    yes/no  (library strategy ncRNA-Seq/miRNA-Seq, or small-RNA wording when strategy empty/OTHER)
  processed_file  first supplementary file that could hold sequences+counts (txt/tsv/csv/fa/fasta), else ''
  genotype_class  wild_type | mutant_or_transgenic | other_line (NIL, breeding line) ; with matched term
  organ           normalized organ/tissue category
  treatment_class control_or_none | treated
Rules are explicit regular expressions so decisions are reviewable.
"""
import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from species_profiles import PROFILES  # noqa: E402

# (?<![a-z]) / (?![a-z]) instead of \b so that "_" and "os" prefixes still match (e.g. KD_CLSY3, OE_Rep1, osnrpd1ab)
MUTANT = re.compile(r'(mutant|mutants|mutation|knock[- ]?out|knock[- ]?down|(?<![a-z])(ko|kd|oe|ox)(?![a-z])|(?<![a-z])rnai(?![a-z])|'
                    r'transgenic|over[- ]?express\w*|amir|crispr|cas9|talen|t-dna|null|insertion|'
                    r'dcl\d|ago\d|rdr\d|nrpd|nrpe|hen1|hst\d?|msp1|eat1|tdl1|pair\d|osdcl|osago|osrdr|clsy|ddm1|'
                    r'drm\d|cmt\d|met1|kyp|suvh\d|ros1|osdrm|fem\d|des1|'
                    r'::|complementation|complemented|silenc\w*|antisense|ami-?r|edited|copies|\d+ ?copy)', re.I)
WILD = re.compile(r'(wild[- ]?type|(?<![a-z])wt(?![a-z])|control plants?|non-?transgenic)', re.I)
OTHER_LINE = re.compile(r'\b(nil|near[- ]isogenic|introgression|ril|cms|male[- ]sterile|maintainer|restorer|hybrid|'
                        r'f1\b|\d+s\b|ptgms|tgms|egms)', re.I)
SRNA_WORDS = re.compile(r'small ?rna|srna|mirna|micro ?rna|sirna|phasirna|ncrna', re.I)
STRESS = re.compile(r'infect|xoo|pathogen|blast|virus|viroid|rsv|botrytis|pseudomonas|phytophthora|fusarium|nematode|'
                    r'meloidogyne|cladosporium|xanthomonas|ralstonia|tylcv|tomv|tmv|cmv|pstvd|drought|cold|salt|heat|chill|freez|stress|'
                    r'hpi\b|\d+ ?h(our)?s? (post|after)|\bn0\b|low n|deficien|zn[- ]|arsenic|as\(iii\)|submerg|anoxi|hypoxi|'
                    r'glyphosate|herbicide|aza|drug|uv|light exposure|dark|wound|insect|planthopper|treated(?! with water)', re.I)
CONTROL = re.compile(r'\b(control|untreated|no treatment|mock|water|none|normal|ck\b|0 ?h\b|non-?treated)', re.I)
def organ_of(r, organs):
    """Organ from the explicit tissue/source fields first; title and stage text only as fallback."""
    parts = r['characteristics'].split(' | ')
    explicit = ' '.join(p.split(':', 1)[1] for p in parts if re.match(r'\s*(tissue|organ|cell type|source_name)\s*:', p, re.I))
    for text in (explicit, r['source'], r['title'],
                 ' '.join(p for p in parts if re.match(r'\s*(developmental stage|stage|age)\s*:', p, re.I))):
        text = text.lower()
        if not text.strip():
            continue
        for name, pat in organs:
            if name == 'mixed_tissues' and text is not explicit.lower():
                continue  # only call "mixed" from an explicit multi-tissue field
            if re.search(pat, text):
                return name
    return 'unclassified'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inventory_dir', nargs='?', default='inventory/geo')
    ap.add_argument('--species', default='rice', choices=sorted(PROFILES))
    a = ap.parse_args()
    prof = PROFILES[a.species]
    IN = os.path.join(a.inventory_dir, 'geo_samples.tsv')
    OUT = os.path.join(a.inventory_dir, 'geo_samples.annotated.tsv')
    rows = list(csv.DictReader(open(IN), delimiter='\t'))
    out_fields = list(rows[0].keys()) + ['subspecies', 'srna_library', 'processed_file', 'genotype_class',
                                         'genotype_evidence', 'organ', 'treatment_class']
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            text = ' '.join([r['title'], r['source'], r['characteristics']]).lower()
            r['subspecies'] = prof['genotype'](r, text)
            strat = r['library_strategy']
            r['srna_library'] = 'yes' if strat in ('ncRNA-Seq', 'miRNA-Seq') or (
                strat in ('', 'OTHER') and SRNA_WORDS.search(text)) else 'no'
            files = [f for f in r['supplementary_files'].split(';') if f and f != 'NONE']
            cand = [f for f in files if re.search(r'\.(txt|tsv|csv|fa|fasta|fna|tag|seq|counts?)(\.gz|\.bz2)?$', f, re.I)]
            r['processed_file'] = cand[0] if cand else ''
            geno_text = ' '.join(p for p in r['characteristics'].split(' | ')
                                 if re.match(r'(genotype|strain|cultivar|line|variety|ecotype|background)', p, re.I))
            geno_text = (geno_text + ' ' + r['title'] + ' ' + r['source']).strip()
            m = MUTANT.search(geno_text) or (prof['extra_mutant'].search(geno_text) if prof['extra_mutant'] else None)
            if m and not WILD.search(geno_text.split(m.group(0))[0][-40:] + m.group(0)):
                r['genotype_class'], r['genotype_evidence'] = 'mutant_or_transgenic', m.group(0)
            elif OTHER_LINE.search(geno_text) or (prof['extra_other_line'] and prof['extra_other_line'].search(geno_text)):
                hit = OTHER_LINE.search(geno_text) or prof['extra_other_line'].search(geno_text)
                r['genotype_class'], r['genotype_evidence'] = 'other_line', hit.group(0)
            else:
                r['genotype_class'], r['genotype_evidence'] = 'wild_type', (WILD.search(geno_text).group(0) if WILD.search(geno_text) else '')
            r['organ'] = organ_of(r, prof['organs'])
            treat = ' '.join(p for p in r['characteristics'].split(' | ')
                             if re.match(r'\s*(treatment|stress|condition|agent|infection|time)', p, re.I))
            treat_all = ' '.join([treat, r['title']]).lower()
            if STRESS.search(treat_all) and not CONTROL.search(treat):
                r['treatment_class'] = 'treated'
            elif treat and not CONTROL.search(treat) and not STRESS.search(treat_all):
                r['treatment_class'] = 'treated'
            else:
                r['treatment_class'] = 'control_or_none'
            w.writerow(r)
    print('wrote', OUT, file=sys.stderr)


if __name__ == '__main__':
    main()
