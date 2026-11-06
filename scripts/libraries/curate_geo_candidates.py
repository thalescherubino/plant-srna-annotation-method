#!/usr/bin/env python3
"""Annotate rice GEO samples for small-RNA library selection (no downloads of data files).

Adds derived columns to inventory/geo/geo_samples.tsv:
  subspecies      japonica_nipponbare | japonica_other | indica | unknown_sativa | non_target
  srna_library    yes/no  (library strategy ncRNA-Seq/miRNA-Seq, or small-RNA wording when strategy empty/OTHER)
  processed_file  first supplementary file that could hold sequences+counts (txt/tsv/csv/fa/fasta), else ''
  genotype_class  wild_type | mutant_or_transgenic | other_line (NIL, breeding line) ; with matched term
  organ           normalized organ/tissue category
  treatment_class control_or_none | treated
Rules are explicit regular expressions so decisions are reviewable.
"""
import csv
import re
import sys

IN = (sys.argv[1] if len(sys.argv) > 1 else 'inventory/geo') + '/geo_samples.tsv'
OUT = (sys.argv[1] if len(sys.argv) > 1 else 'inventory/geo') + '/geo_samples.annotated.tsv'

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
STRESS = re.compile(r'infect|xoo|pathogen|blast|virus|rsv|drought|cold|salt|heat|chill|freez|stress|'
                    r'hpi\b|\d+ ?h(our)?s? (post|after)|\bn0\b|low n|deficien|zn[- ]|arsenic|as\(iii\)|submerg|anoxi|hypoxi|'
                    r'glyphosate|herbicide|aza|drug|uv|light exposure|dark|wound|insect|planthopper|treated(?! with water)', re.I)
CONTROL = re.compile(r'\b(control|untreated|no treatment|mock|water|none|normal|ck\b|0 ?h\b|non-?treated)', re.I)
ORGANS = [
    ('mixed_tissues', r',.*\b(and|,)\b.*|inflorescences, leaves|leaves, roots|mixed|pooled tissues'),
    ('anther_pollen_meiocyte', r'anther|pollen|meiocyte|microspore|tetrad|stamen'),
    ('root', r'\broots?\b|radicle|root tip'),
    ('panicle_inflorescence', r'panicle|inflorescence|spikelet|floret|flower|young ear|rachis|glume|pistil|ovar'),
    ('seed_embryo_endosperm', r'embryo|endosperm|\bseeds?\b|grain|caryopsis|\bdap\b|aleurone'),
    ('callus_culture', r'callus|cell culture|suspension|protoplast'),
    ('leaf', r'flag ?leaf|leaf|leaves|blade'),
    ('sheath_stem_node', r'sheath|stem|culm|node|internode|tiller|crown|shoot base|shoot apex|\bsam\b|meristem'),
    ('seedling_shoot', r'seedlings?|shoot|aerial|whole plant|plantlet|coleoptile'),
]


def organ_of(r):
    """Organ from the explicit tissue/source fields first; title and stage text only as fallback."""
    parts = r['characteristics'].split(' | ')
    explicit = ' '.join(p.split(':', 1)[1] for p in parts if re.match(r'\s*(tissue|organ|cell type)\s*:', p, re.I))
    for text in (explicit, r['source'], r['title'],
                 ' '.join(p for p in parts if re.match(r'\s*(developmental stage|stage|age)\s*:', p, re.I))):
        text = text.lower()
        if not text.strip():
            continue
        for name, pat in ORGANS:
            if name == 'mixed_tissues' and text is not explicit.lower():
                continue  # only call "mixed" from an explicit multi-tissue field
            if re.search(pat, text):
                return name
    return 'unclassified'


def subspecies(r, text):
    org = r['organism'].lower()
    if 'oryza sativa' not in org:
        return 'non_target'
    if re.search(r'nipp?on?bare|nipponbare|nippnbare', text):
        return 'japonica_nipponbare'
    if 'indica' in org or re.search(r'indica|93-?11|9311|ir24|ir64|irbb|pusa basmati|pb-?1\b|zhenshan|minghui|teqing|kasalath|co39|n22|aus\b', text):
        return 'indica_or_aus'
    if 'japonica' in org or re.search(r'japonica|japanica|japponica|kitaake|zhonghua|dongjin|hwayoung|nongken|taichung|'
                                      r'koshihikari|wuyunjing|ilmi|illmi|jowoon|gopum|azucena|\baz_|9522|zhonghua ?11|'
                                      r'heukgyeong|nipponbare', text):
        return 'japonica_other'
    return 'unknown_sativa'


def main():
    rows = list(csv.DictReader(open(IN), delimiter='\t'))
    out_fields = list(rows[0].keys()) + ['subspecies', 'srna_library', 'processed_file', 'genotype_class',
                                         'genotype_evidence', 'organ', 'treatment_class']
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            text = ' '.join([r['title'], r['source'], r['characteristics']]).lower()
            r['subspecies'] = subspecies(r, text)
            strat = r['library_strategy']
            r['srna_library'] = 'yes' if strat in ('ncRNA-Seq', 'miRNA-Seq') or (
                strat in ('', 'OTHER') and SRNA_WORDS.search(text)) else 'no'
            files = [f for f in r['supplementary_files'].split(';') if f and f != 'NONE']
            cand = [f for f in files if re.search(r'\.(txt|tsv|csv|fa|fasta|fna|tag|seq|counts?)(\.gz|\.bz2)?$', f, re.I)]
            r['processed_file'] = cand[0] if cand else ''
            geno_text = ' '.join(p for p in r['characteristics'].split(' | ')
                                 if re.match(r'(genotype|strain|cultivar|line|variety|ecotype|background)', p, re.I))
            geno_text = (geno_text + ' ' + r['title'] + ' ' + r['source']).strip()
            m = MUTANT.search(geno_text)
            if m and not WILD.search(geno_text.split(m.group(0))[0][-40:] + m.group(0)):
                r['genotype_class'], r['genotype_evidence'] = 'mutant_or_transgenic', m.group(0)
            elif OTHER_LINE.search(geno_text):
                r['genotype_class'], r['genotype_evidence'] = 'other_line', OTHER_LINE.search(geno_text).group(0)
            else:
                r['genotype_class'], r['genotype_evidence'] = 'wild_type', (WILD.search(geno_text).group(0) if WILD.search(geno_text) else '')
            r['organ'] = organ_of(r)
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
