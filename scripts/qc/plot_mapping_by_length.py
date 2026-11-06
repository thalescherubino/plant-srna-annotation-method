#!/usr/bin/env python3
"""Plot exact-match (0-mismatch) genome mapping rate by small-RNA length and cultivar group.

Input : exact_mapping_qc_by_length.tsv (one row per library; read-weighted fraction per length 18-30 nt)
Output: <prefix>.png, <prefix>.svg and <prefix>.summary.tsv (median and IQR per group x length)
Palette: reference categorical slots 1-4, validated (light, adjacent pairs): CVD dE >= 9.1, normal-vision dE >= 22.9;
slots 3-4 are below 3:1 contrast, so every line is directly labelled and a summary table is written.
"""
import csv
import statistics
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

inp, prefix = sys.argv[1], sys.argv[2]
REF_LABEL = sys.argv[3] if len(sys.argv) > 3 else 'reference genome'  # e.g. 'Rice IRGSP-1.0 genome'
GROUPS = [  # (key, label, color) in fixed categorical order
    ('japonica_nipponbare', 'Nipponbare (reference cultivar)', '#2a78d6'),
    ('japonica_other', 'Other japonica', '#eb6834'),
    ('indica_or_aus', 'Indica / aus', '#1baf7a'),
    ('unknown_sativa', 'Cultivar not recorded', '#eda100'),
]
INK, INK2, MUTED, GRID, SURFACE = '#0b0b0b', '#52514e', '#8a8984', '#e7e6e2', '#fcfcfb'
LENGTHS = list(range(18, 25))  # plotted: 18-24 nt (user request 2026-11-06); axis still extends to 25 nt
AXIS_TICKS = list(range(18, 26))
MIN_READS = 1000  # a library contributes to a length only if it has >= this many reads of that length

rows = list(csv.DictReader(open(inp), delimiter='\t'))
summary = []
fig, ax = plt.subplots(figsize=(11, 7.5), dpi=200)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
for n in (21, 22, 24):
    ax.axvline(n, color=GRID, lw=1.2, zorder=0)
    ax.text(n, 100.4, f'{n} nt', ha='center', va='bottom', fontsize=12, color=INK2)
for key, label, color in GROUPS:
    libs = [r for r in rows if r['subspecies'] == key]
    if not libs:
        continue
    xs, med, q1, q3 = [], [], [], []
    for L in LENGTHS:
        vals = [100 * float(r[f'frac_mapped_{L}nt']) for r in libs
                if r[f'frac_mapped_{L}nt'] not in ('', None) and int(r[f'reads_{L}']) >= MIN_READS]
        if len(vals) < 2:
            continue
        qs = statistics.quantiles(vals, n=4, method='inclusive')
        xs.append(L); med.append(statistics.median(vals)); q1.append(qs[0]); q3.append(qs[2])
        summary.append({'group': label, 'length_nt': L, 'n_libraries': len(vals), 'median_pct_mapped': round(med[-1], 2),
                        'q1_pct': round(qs[0], 2), 'q3_pct': round(qs[2], 2)})
    ax.fill_between(xs, q1, q3, color=color, alpha=0.10, lw=0, zorder=1)
    ax.plot(xs, med, color=color, lw=2.5, solid_joinstyle='round', solid_capstyle='round', zorder=3,
            label=f'{label} (n={len(libs)})')
    ax.plot(xs, med, 'o', ms=8, color=color, markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=4)

# identity is carried by the legend (no end-of-line labels, per user request); the summary TSV is the table view
ax.set_xlim(17.6, 25.4)
ax.set_ylim(55, 100)
ax.set_xticks(AXIS_TICKS)
ax.set_xlabel('Small-RNA length (nt)', fontsize=14, color=INK2)
ax.set_ylabel('Reads mapping to the reference\nwith 0 mismatches (%)', fontsize=14, color=INK2)
ax.yaxis.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
for s in ('left', 'bottom'):
    ax.spines[s].set_color(MUTED)
ax.tick_params(colors=INK2, labelsize=13)
fig.suptitle('Exact-match mapping rate by small-RNA length and cultivar group', x=0.06, ha='left',
             fontsize=17, color=INK, y=0.985)
ax.set_title(f'Line: median across libraries; band: interquartile range. {REF_LABEL},\n'
             f'Bowtie 1 (-v 0), read-weighted; {len(rows)} libraries; 18–24 nt plotted.', loc='left', fontsize=12,
             color=INK2, pad=26)
leg = ax.legend(loc='lower left', frameon=False, fontsize=13, ncol=2, handlelength=2.2, columnspacing=1.6)
for t in leg.get_texts():
    t.set_color(INK)
fig.tight_layout()
fig.savefig(prefix + '.png', facecolor=SURFACE)
fig.savefig(prefix + '.svg', facecolor=SURFACE)
with open(prefix + '.summary.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(summary)
print('wrote', prefix + '.png', prefix + '.svg', prefix + '.summary.tsv')
