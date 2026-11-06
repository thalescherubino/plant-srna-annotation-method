#!/usr/bin/env python3
"""Select up to N small-RNA libraries for annotation, maximizing organ diversity.

Inputs: annotated GEO samples + format probe (see curate_geo_candidates.py, probe_geo_formats.py).
Rules (2026-11-06 user decisions):
  * rice sRNA libraries only; read-level processed files only (sequence+count, collapsed FASTA, reads, FASTQ);
  * exclude mutants/transgenics and obvious deleterious lines (male-sterile, CMS, TGMS/PTGMS);
    other special lines (NIL, hybrids) are excluded by default and listed;
  * pre-filtered files allowed, flagged;
  * maximize organ diversity: round-robin over organ categories;
  * within an organ: Nipponbare > other japonica > unknown cultivar > indica/aus; control > treated;
    then series with fewer libraries already chosen; then GSM order (deterministic);
  * at most --per-series libraries from one series; total <= --max.
Writes selection, the ranked eligible pool and the exclusion reasons.
"""
import argparse
import collections
import csv
import re

READ_LEVEL = {'tagcount_seq_count', 'tagcount_count_seq', 'fasta_collapsed', 'fasta_reads', 'fastq'}
DELETERIOUS_LINE = re.compile(r'male[- ]sterile|\bcms\b|tgms|ptgms|egms|\b\d+s\b|sterile', re.I)
SUB_RANK = {'japonica_nipponbare': 0, 'japonica_other': 1, 'unknown_sativa': 2, 'indica_or_aus': 3}
NONSTANDARD = re.compile(r'csrna|\[pare\]|(?<![a-z])pare(?![a-z])|degradome|trm-srna|ribo-?z|ribozero|\bip\b|immunoprecip', re.I)
PREFILTER = re.compile(r'match|mapped|map\b|_map|genome|aligned|clean|filtered|uniq', re.I)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', required=True)
    ap.add_argument('--probe', required=True)
    ap.add_argument('--overrides', help='TSV gsm<TAB>field<TAB>value manual corrections (organ, subspecies, exclude)')
    ap.add_argument('--stats', help='library_stats.tsv from inspect_and_standardize_libraries.py (full-file checks)')
    ap.add_argument('--min-reads', type=float, default=1e6)
    ap.add_argument('--allow-special-lines', action='store_true',
                    help='keep non-deleterious special lines (NIL, introgression, F1) ranked after wild type')
    ap.add_argument('--max', type=int, default=100)
    ap.add_argument('--per-series', type=int, default=12)
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    probe = collections.defaultdict(list)
    for p in csv.DictReader(open(a.probe), delimiter='\t'):
        probe[p['gsm']].append(p)
    overrides = collections.defaultdict(dict)
    if a.overrides:
        for o in csv.DictReader(open(a.overrides), delimiter='\t'):
            overrides[o['gsm']][o['field']] = o['value']

    stats = {s['gsm']: s for s in csv.DictReader(open(a.stats), delimiter='\t')} if a.stats else {}
    pool, excluded = [], []
    seen_gsm = {}
    for r in csv.DictReader(open(a.samples), delimiter='\t'):
        for k, v in overrides.get(r['gsm'], {}).items():
            r[k] = v
        reason = None
        text = ' '.join([r['title'], r['source'], r['characteristics']])
        if r['subspecies'] == 'non_target':
            reason = 'not_oryza_sativa'
        elif r['srna_library'] != 'yes':
            reason = 'not_small_rna_library'
        elif r.get('exclude'):
            reason = 'manual:' + r['exclude']
        elif NONSTANDARD.search(text + ' ' + r.get('library_selection', '')):
            reason = 'nonstandard_library:' + NONSTANDARD.search(text + ' ' + r.get('library_selection', '')).group(0)
        elif r['gsm'] in seen_gsm:
            reason = 'duplicate_gsm_in_superseries:' + seen_gsm[r['gsm']]
        elif r['genotype_class'] == 'mutant_or_transgenic':
            reason = 'mutant_or_transgenic:' + r['genotype_evidence']
        elif DELETERIOUS_LINE.search(text):
            reason = 'deleterious_line:' + DELETERIOUS_LINE.search(text).group(0)
        elif r['genotype_class'] == 'other_line' and not a.allow_special_lines:
            reason = 'special_line:' + r['genotype_evidence']
        else:
            files = [p for p in probe.get(r['gsm'], []) if p['format'] in READ_LEVEL]
            if not files:
                fmts = ','.join(sorted({p['format'] for p in probe.get(r['gsm'], [])})) or 'no_probed_file'
                reason = 'no_read_level_file:' + fmts
            else:
                best = sorted(files, key=lambda p: list(READ_LEVEL).index(p['format']) if False else
                              ['tagcount_seq_count', 'tagcount_count_seq', 'fasta_collapsed', 'fasta_reads', 'fastq'].index(p['format']))[0]
                r['selected_file'] = best['file']
                r['file_format'] = best['format']
                r['prefiltered_flag'] = 'yes' if PREFILTER.search(best['file']) else 'unknown'
                r['length_range_probe'] = f"{best['min_len']}-{best['max_len']}"
                s = stats.get(r['gsm'])
                if a.stats and s is None:
                    reason = 'not_inspected'
                elif s is not None and s['integer_counts'] != 'True':
                    reason = 'non_integer_counts(normalized_values)'
                elif s is not None and not s['tagcount_file']:
                    reason = 'no_standard_tagcount_written'
                elif s is not None and float(s['kept_reads']) < a.min_reads:
                    reason = f"low_depth:{s['kept_reads']}"
                elif s is not None and (int(s['min_len']) > 20 or int(s['max_len']) < 24):
                    reason = f"length_truncated:{s['min_len']}-{s['max_len']}"
                elif s is not None and float(s['frac_18_30']) < 0.5:
                    reason = f"few_18_30nt_reads:{s['frac_18_30']}"
                else:
                    if s is not None:
                        r['file_format'] = s['layout']
                        r['length_range_probe'] = f"{s['min_len']}-{s['max_len']}"
                        r['kept_reads'] = s['kept_reads']
                        r['tagcount_file'] = s['tagcount_file']
                        r['frac_21'], r['frac_24'] = s['frac_21'], s['frac_24']
                    seen_gsm[r['gsm']] = r['gse']
                    pool.append(r)
                    continue
        r['exclusion_reason'] = reason
        excluded.append(r)

    # identical standardized content deposited under different GSMs: keep the lowest GSM number
    by_content = collections.defaultdict(list)
    for r in pool:
        s = stats.get(r['gsm'])
        if s is not None and s['tagcount_sha256']:
            by_content[s['tagcount_sha256']].append(r)
    for group in by_content.values():
        if len(group) > 1:
            group.sort(key=lambda r: int(r['gsm'][3:]))
            for dup in group[1:]:
                dup['exclusion_reason'] = 'duplicate_content_of:' + group[0]['gsm']
                excluded.append(dup)
                pool.remove(dup)

    def rank(r):
        return (SUB_RANK.get(r['subspecies'], 9), 0 if r['genotype_class'] == 'wild_type' else 1,
                0 if r['treatment_class'] == 'control_or_none' else 1, r['gsm'])

    by_organ = collections.defaultdict(list)
    for r in sorted(pool, key=rank):
        by_organ[r['organ']].append(r)
    organs = sorted(by_organ, key=lambda o: (o == 'unclassified', o))
    chosen, per_series = [], collections.Counter()
    progress = True
    while len(chosen) < a.max and progress:
        progress = False
        for o in organs:
            if len(chosen) >= a.max:
                break
            cands = [r for r in by_organ[o] if per_series[r['gse']] < a.per_series and 'chosen' not in r]
            if not cands:
                continue
            cands.sort(key=lambda r: rank(r)[:3] + (per_series[r['gse']], r['gsm']))
            r = cands[0]
            r['chosen'] = 'yes'
            per_series[r['gse']] += 1
            r['selection_order'] = len(chosen) + 1
            chosen.append(r)
            progress = True

    keep = ['selection_order', 'gse', 'gsm', 'title', 'organ', 'subspecies', 'treatment_class', 'genotype_class',
            'file_format', 'selected_file', 'prefiltered_flag', 'length_range_probe', 'kept_reads', 'frac_21', 'frac_24',
            'tagcount_file', 'source', 'characteristics']
    with open(a.out_prefix + '.selected.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keep, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(chosen)
    with open(a.out_prefix + '.eligible_pool.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keep[1:] + ['chosen'], delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(sorted(pool, key=lambda r: (r['organ'], rank(r))))
    with open(a.out_prefix + '.excluded.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['gse', 'gsm', 'title', 'organ', 'subspecies', 'exclusion_reason'],
                           delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(excluded)

    C = collections.Counter
    print(f'eligible={len(pool)} selected={len(chosen)} excluded={len(excluded)}')
    print('eligible by organ:', dict(C(r['organ'] for r in pool)))
    print('selected by organ:', dict(C(r['organ'] for r in chosen)))
    print('selected by subspecies:', dict(C(r['subspecies'] for r in chosen)))
    print('selected by format:', dict(C(r['file_format'] for r in chosen)))
    print('selected series:', len({r['gse'] for r in chosen}), dict(per_series.most_common(8)))
    print('exclusion reasons:', C(r['exclusion_reason'].split(':')[0] for r in excluded).most_common())


if __name__ == '__main__':
    main()
