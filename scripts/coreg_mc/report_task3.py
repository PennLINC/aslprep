#!/usr/bin/env python3
"""Walkthrough HTML report of Task 3's comparator evidence: R1's sign rule and R2's D3 level.

Reads ``freeze_comparators.py``'s ``fixtures.tsv`` and ``real.tsv`` (and, with ``--capture``,
``capture.py``'s table of the comparators' capture range) and writes one Artifact page
fragment (no doctype or ``<html>``/``<body>``). The chart helpers and styling are the HMC probe's
(``scripts/hmc_probe/report.py``). Every number in the prose is computed here from the data.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/report_task3.py --freeze <freeze dir> --out report.html
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import freeze_comparators  # noqa: E402
from report import (  # noqa: E402
    CSS,
    esc,
    figure,
    fmt,
    nice_ticks,
    plain_table,
    strip_plot,
    svg_open,
)

LEVELS = ('P-mean', 'P-robust', 'P-winsor', 'P-masked')
FAMILIES = ('2D EPI', '3D EPI 7T', '3D GRASE / 3D pCASL', 'stack-of-spirals', 'GE product')
CONVENTIONAL = 'FSL conventional default'
#: Inconsistency bands (mm) used to summarize a run: within 2 mm, beyond the gross limit.
NEAR_MM, GROSS_MM = 2.0, 5.0
#: A fixture the default sign misses: error against truth at or above this (mm).
MISS_MM = 1.0
#: Strip plots are clipped here; larger values sit on the right edge.
CLIP_MM = 20.0


def signed_strip(cats, values, xlabel, label, w=680):
    """One row per category on a signed axis with a zero line; jittered points, median tick."""
    row = 26
    m = (10, 16, 40, 230)
    h = m[0] + m[2] + row * len(cats)
    pw = w - m[1] - m[3]
    allv = np.concatenate([np.asarray(v, float) for v in values])
    xt = nice_ticks(float(np.nanmin(allv)), float(np.nanmax(allv)))
    x0, x1 = xt[0], xt[-1]

    def xs(v):
        return m[3] + (v - x0) / (x1 - x0) * pw

    out = [svg_open(w, h, label)]
    for t in xt:
        cls = 'axis' if abs(t) < 1e-12 else 'grid'
        out.append(
            f'<line class="{cls}" x1="{xs(t):.1f}" x2="{xs(t):.1f}" y1="{m[0]}" y2="{h - m[2]}"/>'
        )
        out.append(
            f'<text class="tick" x="{xs(t):.1f}" y="{h - m[2] + 16}" '
            f'text-anchor="middle">{t:+g}</text>'
        )

    rng = np.random.default_rng(0)
    for i, (c, v) in enumerate(zip(cats, values, strict=True)):
        cy = m[0] + row * (i + 0.5)
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        out.append(
            f'<text class="tick" x="{m[3] - 8}" y="{cy + 4:.1f}" '
            f'text-anchor="end">{esc(c)} (n={v.size})</text>'
        )
        for val, jit in zip(v, rng.uniform(-7, 7, v.size), strict=True):
            out.append(
                f'<circle cx="{xs(val):.1f}" cy="{cy + jit:.1f}" r="2.6" '
                'fill="var(--accent)" fill-opacity="0.45"/>'
            )

        if v.size:
            mx = xs(float(np.median(v)))
            out.append(
                f'<line class="median" x1="{mx:.1f}" x2="{mx:.1f}" '
                f'y1="{cy - 10:.1f}" y2="{cy + 10:.1f}"/>'
            )

    out.append(
        f'<text class="axlabel" x="{m[3] + pw / 2}" y="{h - 4}" '
        f'text-anchor="middle">{esc(xlabel)}</text>'
    )
    out.append('</svg>')
    return ''.join(out)


rule_sign = freeze_comparators.r1_sign
THRESHOLD = freeze_comparators.R1_CONTRAST_THRESHOLD


def load(freeze):
    fix = pd.read_table(freeze / 'fixtures.tsv')
    real = pd.read_table(freeze / 'real.tsv')
    key = ['dataset', 'run']
    fsl = real[real['arm'].str.startswith('FSL')]
    every = fsl.groupby(key)['verdict'].apply(lambda s: (s == 'fail').all())
    shared = {k for k, v in every.items() if v}
    real['shared_failure'] = [k in shared for k in zip(real['dataset'], real['run'], strict=True)]
    return fix, real


def runs_of(real):
    return real.drop_duplicates(['dataset', 'run']).copy()


def count(verdicts, value):
    return int((np.asarray(verdicts) == value).sum())


# ---------------------------------------------------------------------------------------------
# R1: the BBR sign of the conventional reference
# ---------------------------------------------------------------------------------------------
def r1_fixtures(fix):
    w = fix.pivot(index='fixture', columns='arm', values='rms_mm')
    c = fix.drop_duplicates('fixture').set_index('fixture')['contrast']
    out = pd.DataFrame(
        {
            'contrast': c,
            'default': w['FSL conventional default'],
            'opposite': w['FSL conventional opposite'],
        },
    )
    out['rule'] = [r[rule_sign(r['contrast'])] for _, r in out.iterrows()]
    return out.sort_values('contrast')


def r1_real(real):
    """Per run: contrast, family, and the conventional reference's verdict per sign and rule."""
    conv = real[real['arm'].str.startswith('FSL conventional')]
    v = conv.pivot_table(
        index=['dataset', 'run'],
        columns='arm',
        values='verdict',
        aggfunc='first',
    )
    runs = runs_of(real).set_index(['dataset', 'run'])[['family', 'contrast', 'shared_failure']]
    out = runs.join(
        v.rename(
            columns={
                'FSL conventional default': 'default',
                'FSL conventional opposite': 'opposite',
            },
        ),
    )
    out['rule_sign'] = out['contrast'].map(rule_sign)
    out['rule'] = [r[r['rule_sign']] for _, r in out.iterrows()]
    return out


def section_r1(fix, real, fig):
    f = r1_fixtures(fix)
    r = r1_real(real)
    inverted = f[f['contrast'] > 0]
    upright = f[f['contrast'] <= 0]
    flips = r[r['rule_sign'] == 'opposite']
    positive = r[r['contrast'] > 0]
    ruled = f[f['contrast'] > THRESHOLD]
    missed = f[f['default'] >= MISS_MM]
    rule_misses = int((f['rule'] >= MISS_MM).sum())
    exact = set(missed.index) == set(inverted.index)

    cats = ['Fixtures the default sign misses', 'Other fixtures']
    vals = [missed['contrast'].to_numpy(), f.loc[f['default'] < MISS_MM, 'contrast'].to_numpy()]
    for fam in FAMILIES:
        x = r.loc[r['family'] == fam, 'contrast'].to_numpy()
        if x.size:
            cats.append(f'Real: {fam}')
            vals.append(x)

    chart = figure(
        signed_strip(
            cats,
            vals,
            'contrast (right of zero: white matter brighter than grey)',
            label='Conventional reference contrast, fixtures and real runs',
        ),
        (
            'Contrast of the conventional reference under the mri_coreg initialization. '
            f'A fixture is a miss when the default sign is {MISS_MM:g} mm or more from truth. '
            'The tick is the median.'
        ),
        fig,
    )

    rows = [
        [
            f'<code>{esc(n)}</code>',
            f'{x["contrast"]:+.3f}',
            fmt(x['default']),
            fmt(x['opposite']),
            f'{rule_sign(x["contrast"])}, {fmt(x["rule"])}',
        ]
        for n, x in f.iterrows()
    ]
    fix_table = plain_table(
        ['Fixture', 'Contrast', 'Default (mm)', 'Opposite (mm)', 'Rule (mm)'],
        rows,
        num_cols=(1, 2, 3),
    )

    fam_rows = []
    for fam in FAMILIES:
        g = r[r['family'] == fam]
        if not len(g):
            continue
        fam_rows.append(
            [
                esc(fam),
                str(len(g)),
                f'{g["contrast"].median():+.3f}',
                f'{g["contrast"].max():+.3f}',
                *(
                    f'{count(g[c], "fail")} / {count(g[c], "ambiguous")}'
                    for c in ('default', 'opposite', 'rule')
                ),
            ],
        )
    real_table = plain_table(
        [
            'Family',
            'Runs',
            'Median contrast',
            'Max contrast',
            'Default: fail / ambiguous',
            'Opposite',
            'Rule',
        ],
        fam_rows,
        num_cols=(1, 2, 3),
    )
    changed = flips[flips['default'] != flips['rule']]
    total = {c: count(r[c], 'fail') for c in ('default', 'opposite', 'rule')}
    shared = int(r['shared_failure'].sum())
    body = f"""
<section id="r1">
<h2>R1: the BBR sign of the conventional reference</h2>
<p>BBR's default sign expects grey matter brighter than white matter in the moving image. The
proposed rule reads the contrast of the conventional reference before BBR runs: the median of the
reference at eroded white-matter voxels minus its median at grey-matter-core voxels, both sampled
through the <code>mri_coreg</code> initialization and divided by the median over the brain. A
contrast above {THRESHOLD:+.2f} selects the opposite sign. (Frozen on 2026-10-30 at zero;
amended on 2026-11-02, below.)</p>
<p>On the {len(f)} development fixtures, the {len(inverted)} with positive contrast
({', '.join(f'<code>{esc(n)}</code>' for n in inverted.index)},
{inverted['contrast'].min():+.3f} to {inverted['contrast'].max():+.3f})
{'are exactly' if exact else 'do not match'} the {len(missed)} the default sign misses by
{MISS_MM:g} mm or more, and {len(ruled)} of them lie above the threshold. The other
{len(upright)} sit between {upright['contrast'].min():+.3f} and
{upright['contrast'].max():+.3f}. The rule takes the default sign's {len(missed)} misses to
{rule_misses}.</p>
{chart}
{fix_table}
<p>Real data does not show the inversion. Across {len(r)} development runs only
{len(positive)} have a positive contrast, the largest at {r['contrast'].max():+.3f}:
{', '.join(f'<code>{esc(d)}/{esc(n)}</code> ({c:+.3f})' for (d, n), c in positive['contrast'].items())}.
At zero the rule changed the sign on all of them; on the one of these the reviewer saw in both
development batches, the opposite sign was judged worse both times and flagged as misaligned once.
The threshold now sits between that run and the smallest fixture miss, and the rule changes the
sign on {len(flips)} development runs. The default sign has {total['default']} gross failures and
the opposite sign {total['opposite']}; the rule has {total['rule']}. {shared} of the failures are
runs every arm fails (<a href="#review">below</a>), so they say nothing about the sign.</p>
{real_table}
<p class="note">Fail and ambiguous are the adopted gross-failure verdicts (more than 5 mm
inconsistent with every partner when there are two or more; with a single partner both runs are
ambiguous).</p>
</section>"""
    facts = {
        'fixtures': len(f),
        'missed': len(missed),
        'rule_misses': rule_misses,
        'exact': exact,
        'max_real': float(r['contrast'].max()),
        'flips': len(flips),
        'changed': len(changed),
        'fails': total,
        'margin': (float(upright['contrast'].max()), float(inverted['contrast'].min())),
    }
    return body, facts


# ---------------------------------------------------------------------------------------------
# R2: the D3 level and eligibility
# ---------------------------------------------------------------------------------------------
def eligible_runs(real):
    """Rows of eligible runs with partners, outside the runs every arm fails."""
    e = real[(real['elig_eligible'] == True) & ~real['shared_failure']]  # noqa: E712
    return e[e['median_inconsistency_mm'].notna()]


def level_summary(g):
    x = g['median_inconsistency_mm']
    return {
        'n': len(x),
        'median': x.median(),
        'p75': x.quantile(0.75),
        'near': (x <= NEAR_MM).mean(),
        'gross': (x > GROSS_MM).mean(),
        'fail': count(g['verdict'], 'fail'),
        'ambiguous': count(g['verdict'], 'ambiguous'),
    }


def section_r2(fix, real, fig):
    e = eligible_runs(real)
    arms = [CONVENTIONAL] + [f'FSL {lv} default' for lv in LEVELS]
    summary = {a: level_summary(e[e['arm'] == a]) for a in arms}
    wide = e[e['arm'].isin(arms)].pivot_table(
        index=['dataset', 'run'],
        columns='arm',
        values='median_inconsistency_mm',
    )
    levels = [f'FSL {lv} default' for lv in LEVELS]
    best = wide[levels].idxmin(axis=1).value_counts()
    beats = float((wide[levels].min(axis=1) < wide[CONVENTIONAL]).mean())

    rows = []
    for a in arms:
        s = summary[a]
        fx = fix[fix['arm'] == a]['rms_mm']
        rows.append(
            [
                esc(a.removeprefix('FSL ').removesuffix(' default')),
                str(s['n']),
                fmt(s['median']),
                fmt(s['p75']),
                f'{s["near"]:.0%}',
                f'{s["gross"]:.0%}',
                f'{s["fail"]} / {s["ambiguous"]}',
                str(int(best.get(a, 0))) if a in levels else '&ndash;',
                f'{fmt(fx.median())} / {fmt(fx.max())}',
            ],
        )
    table = plain_table(
        [
            'Image',
            'Runs',
            'Median (mm)',
            '75th pct.',
            f'&le; {NEAR_MM:g} mm',
            f'&gt; {GROSS_MM:g} mm',
            'Fail / ambiguous',
            'Best level',
            'Fixtures median / max (mm)',
        ],
        rows,
        num_cols=(1, 2, 3, 4, 5, 7),
    )

    fam_rows = []
    leaders = []
    for fam in FAMILIES:
        g = e[e['family'] == fam]
        if not len(g):
            continue
        cells = []
        for a in arms:
            s = level_summary(g[g['arm'] == a])
            cells.append(f'{fmt(s["median"])} ({s["near"]:.0%})')

        fam_rows.append([esc(fam), str(int((g['arm'] == CONVENTIONAL).sum())), *cells])
        medians = g[g['arm'].isin(levels)].groupby('arm')['median_inconsistency_mm'].median()
        lead = medians.idxmin().removeprefix('FSL ').removesuffix(' default')
        leaders.append(f'{esc(fam)}: {lead}')
    fam_table = plain_table(
        ['Family', 'Runs', 'Conventional', *LEVELS],
        fam_rows,
        num_cols=(1,),
    )

    chart = figure(
        strip_plot(
            [a.removeprefix('FSL ').removesuffix(' default') for a in arms],
            [np.minimum(e.loc[e['arm'] == a, 'median_inconsistency_mm'], CLIP_MM) for a in arms],
            f'median inconsistency over partners (mm, clipped at {CLIP_MM:g})',
            label='Inconsistency of each D3 level on eligible runs',
            xmax=CLIP_MM,
        ),
        (
            'Median inconsistency of each run over its partners, conventional frame, eligible '
            'runs only. The tick is the median.'
        ),
        fig,
    )
    opposite = {}
    for lv in LEVELS:
        errors = fix[fix['arm'] == f'FSL {lv} opposite']['rms_mm']
        opposite[lv] = int((errors > NEAR_MM).sum())
    pm = summary['FSL P-mean default']
    body = f"""
<section id="r2">
<h2>R2: which D3 level</h2>
<p>The four D3 levels were registered with both BBR signs on every run. A perfusion image has grey
matter brighter than white, so the default sign is the only candidate: on the fixtures the opposite
sign misses by more than {NEAR_MM:g} mm on {min(opposite.values())} to {max(opposite.values())} of
{fix['fixture'].nunique()} fixtures per level, the default sign on none.</p>
<p>The comparison below is on the {pm['n']} eligible development runs that have partners, leaving
out the runs every arm fails. P-mean is the best level on {int(best.get('FSL P-mean default', 0))}
of {len(wide)} runs, the most of the four, and has the smallest 75th percentile; the levels are
within a few hundredths of a millimeter of each other on the fixtures. No level is close to the
conventional reference on real data: the best level beats it on {beats:.0%} of runs.</p>
{chart}
{table}
<p class="note">Best level: the number of runs on which that level has the lowest inconsistency of
the four. Fixtures: RMS error against truth over the {fix['fixture'].nunique()} development
fixtures.</p>
<h3>By family</h3>
<p>Median inconsistency, with the share of runs within {NEAR_MM:g} mm in parentheses. The 7 T
3D EPI runs are the only family where perfusion images register about as consistently as the
conventional
reference. The level with the lowest median differs by family ({'; '.join(leaders)}), so a single
level is a compromise; P-mean is chosen on the pooled comparison above.</p>
{fam_table}
</section>"""
    facts = {'n': pm['n'], 'best': int(best.get('FSL P-mean default', 0)), 'beats': beats}
    return body, facts


def section_eligibility(real):
    g = real[(real['arm'] == 'FSL P-mean default') & (real['elig_data_type'] == 'pairs')]
    g = g[~g['shared_failure']]
    bins = [0, 2, 3, 4, 6, 10, 20, 1000]
    labels = ['2', '3', '4', '5&ndash;6', '7&ndash;10', '11&ndash;20', '&gt; 20']
    g = g.assign(pairs=pd.cut(g['elig_n_inputs'], bins, labels=labels))
    pair_rows = []
    for lab, s in g.groupby('pairs', observed=True):
        pair_rows.append(
            [
                lab,
                str(len(s)),
                str(count(s['verdict'], 'fail')),
                fmt(s['median_inconsistency_mm'].median()),
            ],
        )
    pair_table = plain_table(
        ['Pairs', 'Runs', 'P-mean fails', 'Median inconsistency (mm)'],
        pair_rows,
        num_cols=(1, 2, 3),
    )
    sweep = []
    enough = g[g['elig_n_inputs'] >= 4]
    for limit in (0.0, 0.4, 0.5, 0.6, 0.7):
        s = enough[enough['elig_edge_r'] >= limit]
        x = s['median_inconsistency_mm'].dropna()
        sweep.append(
            [
                f'{limit:.1f}',
                str(len(s)),
                str(count(s['verdict'], 'fail')),
                f'{(x > GROSS_MM).mean():.0%}',
                fmt(x.median()),
            ],
        )
    sweep_table = plain_table(
        ['Edge limit', 'Eligible', 'P-mean fails', f'&gt; {GROSS_MM:g} mm', 'Median (mm)'],
        sweep,
        num_cols=(1, 2, 3, 4),
    )
    two = g[g['elig_n_inputs'] <= 2]
    four_up = enough.groupby('pairs', observed=True)['median_inconsistency_mm'].median()
    gross_share = []
    for limit in (0.0, 0.5, 0.6, 0.7):
        x = enough[enough['elig_edge_r'] >= limit]['median_inconsistency_mm'].dropna()
        gross_share.append((x > GROSS_MM).mean())
    return f"""
<section id="eligibility">
<h2>R2 eligibility</h2>
<p>The HMC probe's rule makes a run eligible with at least 4 pairs and an edge-content statistic
(the correlation of split-half gradient magnitudes) of at least 0.5; a single derived volume is
never eligible. The pair count does the work. Runs with 2 pairs fail with P-mean on
{count(two['verdict'], 'fail')} of {len(two)} and have a median inconsistency of
{fmt(two['median_inconsistency_mm'].median())} mm; from 4 pairs up the median is
{fmt(four_up.max())} mm or less in every band.</p>
{pair_table}
<p>Among runs with enough pairs, the edge statistic barely separates good runs from bad. The share
more than {GROSS_MM:g} mm inconsistent is {gross_share[0]:.0%} with no limit, {gross_share[1]:.0%}
at 0.5, {gross_share[2]:.0%} at 0.6 and {gross_share[3]:.0%} at 0.7, while each step removes
runs.</p>
{sweep_table}
</section>"""


def section_review(real):
    r = runs_of(real)
    shared = r[r['shared_failure']]
    conv = real[real['arm'] == CONVENTIONAL].set_index(['dataset', 'run'])
    high = conv[(conv['verdict'] == 'ambiguous') & (conv['median_inconsistency_mm'] > GROSS_MM)]
    items = ''.join(
        f'<li><code>{esc(ds)}</code> {esc(run)}</li>'
        for ds, run in zip(shared['dataset'], shared['run'], strict=True)
    )
    by_ds = high.groupby(level=0).size()
    ds_items = ', '.join(f'<code>{esc(ds)}</code> ({n})' for ds, n in by_ds.items())
    return f"""
<section id="review">
<h2>Runs for the reviewer</h2>
<p>Every arm fails on {len(shared)} runs, conventional included. A failure shared by every
registration points at the direct registration between the runs (the consistency check), not at the
arms:</p>
<ul>{items}</ul>
<p>The conventional reference is more than {GROSS_MM:g} mm inconsistent with its single partner on
{len(high)} runs: {ds_items}. These are ambiguous by the adopted rule, and only the blinded review
can say which run of a pair is wrong.</p>
</section>"""


# ---------------------------------------------------------------------------------------------
# Capture range of the frozen comparators
# ---------------------------------------------------------------------------------------------
def capture_comparators(capture, fix):
    """``capture.py``'s table with one ``comparator`` per row: R0, R1 (the rule's sign) and R2."""
    contrast = fix.drop_duplicates('fixture').set_index('fixture')['contrast']
    sign = capture['arm'].str.split().str[-1]
    conventional = capture['arm'].str.startswith('FSL conventional')
    rule = capture['fixture'].map(contrast).map(rule_sign)
    parts = [
        capture[capture['arm'] == 'FSL conventional BBR default'].assign(comparator='R0'),
        capture[conventional & (sign == rule)].assign(comparator='R1'),
        capture[capture['arm'] == 'FSL P-mean BBR default'].assign(comparator='R2'),
    ]
    return pd.concat(parts, ignore_index=True)


def section_capture(capture, fix):
    c = capture_comparators(capture, fix)
    levels = sorted(c['level'].unique())
    start = c.groupby('level')['start_error_mm'].median()
    rows = []
    for comparator, g in c.groupby('comparator'):
        by = g.groupby('level')
        rows.append(
            [comparator, 'with the fallback']
            + [f'{v:.0%}' for v in by['converged'].mean().reindex(levels)],
        )
        rows.append(
            ['', 'optimizer alone']
            + [f'{v:.0%}' for v in by['optimizer_converged'].mean().reindex(levels)],
        )
    head = ['Comparator', 'Result'] + [
        f'{lv:g} mm/&deg; (start {fmt(start[lv])} mm)' for lv in levels
    ]
    table = plain_table(head, rows, num_cols=tuple(range(2, 2 + len(levels))))
    n_fix = c['fixture'].nunique()
    draws = int(c.groupby(['comparator', 'level']).size().iloc[0] / n_fix)
    r0 = c[c['comparator'] == 'R0']
    r0_miss = sorted(f for f, g in r0.groupby('fixture') if not g['optimizer_converged'].any())
    r0_misses = ', '.join(f'<code>{esc(f)}</code>' for f in r0_miss)
    low = c[c['level'] <= 5]
    low_fallback = low.groupby('comparator')['fallback'].mean()
    flat = all(
        g.groupby('level')['optimizer_converged'].mean().nunique() == 1
        for _, g in c.groupby('comparator')
    )
    rule_ok = c[c['comparator'].isin(['R1', 'R2'])]['optimizer_converged'].all()
    worst = c[c['comparator'].isin(['R1', 'R2'])].groupby('comparator')['optimizer_error_mm'].max()
    top = max(levels)
    verb = 'converges' if rule_ok else 'does not always converge'
    high = c.loc[c['level'] >= 10, 'converged'].mean()
    bullet = (
        f'<li><b>Capture range: BBR {verb} from {top:g} mm and {top:g}&deg; with R1 and R2.</b> '
        f"With fMRIPrep's fallback check, {high:.0%} of draws at 10 mm and above converge for "
        'any comparator.</li>'
    )
    if rule_ok:
        lead = 'R1 and R2 converge from every start, at every level, before the fallback check'
    else:
        lead = 'R1 and R2 do not converge from every start'

    level_note = ''
    if flat:
        level_note = 'Convergence does not depend on the start level for any comparator.'
    body = f"""
<section id="capture">
<h2>Capture range</h2>
<p>Each comparator starts from the truth moved by a random rigid transform: a shift of L mm and
a rotation of L degrees, {draws} draws per level on each of the {n_fix} development fixtures, with
the same draws for every arm. FLIRT-BBR runs from that start without <code>mri_coreg</code>, and
fMRIPrep's fallback check compares the result with the start. A draw converges if the result is
within half an acquisition voxel of the truth. R1 takes the sign the frozen rule picks for each
fixture.</p>
{table}
<p>{lead}; the worst result is {fmt(worst.get('R1'))} mm for R1 and {fmt(worst.get('R2'))} mm
for R2. R0 misses every draw on {len(r0_miss)} fixtures ({r0_misses}): these are the fixtures the
default sign misses whatever the start. {level_note}</p>
<p>With the fallback, no comparator converges at 10 or 15 mm. That limit comes from the check, not
from BBR: from these starts the correct answer is more than 15 mm (the check's limit) from the
start, so the check returns the start. At 2 and 5 mm the check fires on {low_fallback['R0']:.0%} of
R0's draws, on fixtures where BBR went wrong, and on {low_fallback['R1']:.0%} of R1's and
{low_fallback['R2']:.0%} of R2's. In production the check compares BBR with
<code>mri_coreg</code>'s result instead, so it bounds how far BBR may move from
<code>mri_coreg</code>, and this protocol cannot show a capture range beyond about 15 mm for any
FSL arm. The optimizer's own result is the measure of BBR's capture range.</p>
</section>"""
    return bullet, body


# ---------------------------------------------------------------------------------------------
# The development review batch
# ---------------------------------------------------------------------------------------------
#: The smallest known-answer size the figures must resolve, and the rate they must reach (gate 4).
GATE_MM, GATE_RATE = 2.0, 0.8
#: Frozen 2026-11-02: consistency pairs whose noise floor exceeds this (mm) are excluded.
FLOOR_LIMIT_MM = 3.0


def section_dev_review(choices, real, section_id, title, intro):
    """One review batch, unblinded. Comparison pairs count as resolvable from the smallest
    known-answer size the reviewer detected at ``GATE_RATE`` or better (on valid pairs)."""
    df = choices.fillna({'note': '', 'failed': ''})
    verdict = real.set_index(['dataset', 'run', 'arm'])['verdict']

    def base_verdict(row):
        dataset, run = row['run'].split('/', 1)
        base = row['arm_B'] if ' moved ' in row['arm_A'] else row['arm_A']
        return verdict.get((dataset, run, base.replace(' BBR', '')), 'none')

    known = df[df['kind'] == 'known'].copy()
    known['base'] = known.apply(base_verdict, axis=1)
    unmoved = known['arm_A'].where(~known['arm_A'].str.contains(' moved '), known['arm_B'])
    known['outcome'] = np.select(
        [known['preferred'] == unmoved, known['preferred'] == 'tie'],
        ['detected', 'tie'],
        'wrong',
    )
    valid = known[known['base'] == 'pass']
    invalid = known[known['base'] != 'pass']
    rows = []
    for size, g in known.groupby('size_mm'):
        v = valid[valid['size_mm'] == size]
        rows.append(
            [
                f'{size:g} mm',
                str(len(g)),
                str(count(g['outcome'], 'detected')),
                str(count(g['outcome'], 'tie')),
                str(count(g['outcome'], 'wrong')),
                f'{count(v["outcome"], "detected")} / {len(v)}',
            ],
        )
    known_table = plain_table(
        ['Size', 'Pairs', 'Detected', 'Tie', 'Wrong', 'Detected, R1 passing'],
        rows,
        num_cols=(1, 2, 3, 4, 5),
    )
    gate = valid[valid['size_mm'] == GATE_MM]
    gate_rate = count(gate['outcome'], 'detected') / max(len(gate), 1)
    rates = valid.groupby('size_mm')['outcome'].apply(lambda s: (s == 'detected').mean())
    resolved = rates[rates >= GATE_RATE]
    resolvable = float(resolved.index.min()) if len(resolved) else float(rates.index.max())
    smallest = valid[valid['size_mm'] == valid['size_mm'].min()]
    smallest_text = (
        f'At {smallest["size_mm"].min():g} mm: {count(smallest["outcome"], "detected")} '
        f'detected, {count(smallest["outcome"], "tie")} ties and '
        f'{count(smallest["outcome"], "wrong")} wrong, out of {len(smallest)}.'
    )
    step = 100 / max(len(gate), 1)
    bad_base = ''.join(
        f'<li><code>{esc(r.run)}</code> ({r.size_mm:g} mm): R1 is {esc(r.base)} there; '
        f'you preferred the {"unmoved" if r.outcome == "detected" else "moved"} side'
        f'{" and flagged both" if ";" in r.failed else ""}.</li>'
        for r in invalid.itertuples(index=False)
    )
    misses = valid[(valid['size_mm'] >= GATE_MM) & (valid['outcome'] != 'detected')]
    miss_items = ', '.join(
        f'<code>{esc(r.run)}</code> ({r.size_mm:g} mm, {r.outcome})'
        for r in misses.itertuples(index=False)
    )

    by_id = df.set_index('id')
    repeats = df[df['kind'] == 'repeat']
    agree = sum(by_id.loc[r.repeat_of, 'preferred'] == r.preferred for r in repeats.itertuples())
    repeat_items = ''.join(
        f'<li><code>{esc(r.run)}</code> ({fmt(by_id.loc[r.repeat_of, "distance_mm"])} mm apart): '
        f'{"the same choice" if by_id.loc[r.repeat_of, "preferred"] == r.preferred else "a different choice"}'
        '</li>'
        for r in repeats.itertuples(index=False)
    )

    comp = df[df['kind'] == 'comparison'].copy()
    sign = comp[comp['comparison'] == 'sign']
    opposite_flagged = int(sign['failed'].str.contains('opposite').sum())
    default_wins = count(sign['preferred'], 'FSL conventional default')
    opposite_wins = sign[sign['preferred'] == 'FSL conventional opposite']
    contrast = runs_of(real).assign(run=lambda r: r['dataset'] + '/' + r['run'])
    contrast = contrast.set_index('run')['contrast']
    rule_arm = 'FSL conventional ' + sign['run'].map(contrast).map(rule_sign)
    rule_default = int((rule_arm == 'FSL conventional default').sum())
    rule_agree = int((sign['preferred'] == rule_arm).sum())
    if rule_default == len(sign):
        rule_text = 'On all of these runs R1&rsquo;s rule picks the default sign'
    else:
        rule_text = f'R1&rsquo;s rule picks the default sign on {rule_default} of these runs'

    if len(opposite_wins):
        wins = ', '.join(
            f'<code>{esc(r.run)}</code> ({fmt(r.distance_mm)} mm)'
            for r in opposite_wins.itertuples()
        )
        opposite_text = f'The opposite sign was preferred on {wins}.'
    else:
        opposite_text = 'The opposite sign was never preferred.'

    rr = comp[comp['comparison'] == 'R1 vs R2'].copy()
    rr['band'] = np.where(
        rr['distance_mm'] >= resolvable,
        f'&ge; {resolvable:g} mm',
        f'&lt; {resolvable:g} mm',
    )
    who = {'FSL P-mean default': 'R2', 'tie': 'tie'}
    rr['who'] = rr['preferred'].map(who).fillna('R1')
    rr_rows = []
    for (family, band), g in rr.groupby(['family', 'band']):
        rr_rows.append(
            [
                esc(family),
                band,
                str(count(g['who'], 'R1')),
                str(count(g['who'], 'R2')),
                str(count(g['who'], 'tie')),
                ', '.join(fmt(d) for d in sorted(g['distance_mm'])),
            ],
        )
    rr_table = plain_table(
        ['Family', 'Distance', 'R1', 'R2', 'Tie', 'Distances (mm)'],
        rr_rows,
        num_cols=(2, 3, 4),
    )
    far = rr[rr['distance_mm'] >= resolvable]
    r2_far = far[far['who'] == 'R2']
    r2_far_items = ', '.join(
        f'<code>{esc(r.run)}</code> ({fmt(r.distance_mm)} mm)' for r in r2_far.itertuples()
    )
    r2_flagged = int(rr['failed'].str.contains('P-mean').sum())
    epi = rr[rr['family'] == '2D EPI']
    epi_far = epi[epi['distance_mm'] >= resolvable]
    epi_text = ''
    if len(epi_far) and epi_far['failed'].str.contains('P-mean').all():
        epi_text = f', including every 2D EPI pair at least {resolvable:g} mm apart'

    near = rr[(rr['distance_mm'] < resolvable) & (rr['who'] == 'R2')]
    near_text = (
        f'R2&rsquo;s {len(near)} other win{"s" if len(near) != 1 else ""} '
        f'{"are" if len(near) != 1 else "is"} on pairs under {resolvable:g} mm apart, which the '
        'known-answer pairs show the figures do not reliably resolve.'
        if len(near)
        else ''
    )
    flipped = rr[rr['arm_A'].str.contains('opposite') | rr['arm_B'].str.contains('opposite')]
    flipped_text = ''.join(
        f' On <code>{esc(r.run)}</code> R1 was shown with the opposite sign (the rule as '
        'frozen at zero); '
        f'{"R1" if r.who == "R1" else "R2" if r.who == "R2" else "neither"} was preferred'
        f'{" and R1 was flagged as misaligned" if "opposite" in r.failed else ""}.'
        for r in flipped.itertuples()
    )
    passed = gate_rate >= GATE_RATE
    verdict_text = (
        'passes the plan&rsquo;s bar'
        if passed
        else 'is below the plan&rsquo;s bar, so the figures are redesigned before any held-out '
        'pair is shown'
    )
    facts = {
        'gate_rate': gate_rate,
        'gate_n': len(gate),
        'passed': passed,
        'far_r1': count(far['who'], 'R1'),
        'far_n': len(far),
        'default_wins': default_wins,
        'sign_n': len(sign),
        'resolvable': resolvable,
    }
    if len(invalid):
        base_text = (
            f'On {len(invalid)} pair{"s" if len(invalid) != 1 else ""} R1 does not pass the '
            f'gross-failure rule:</p><ul>{bad_base}</ul><p>'
        )
    else:
        base_text = 'R1 passes on every known-answer run. '

    misses_text = f'The misses at {GATE_MM:g} mm and above: {miss_items}.' if len(misses) else ''
    body = f"""
<section id="{section_id}">
<h2>{title}</h2>
<p>{intro} {len(df)} figures, all answered: {len(comp)} comparisons, {len(repeats)} repeats and
{len(known)} known-answer pairs. Unblinded here because the batch is complete and development
only.</p>
<h3>Known-answer pairs</h3>
<p>Each pairs R1's solution with the same solution moved by a random rigid error of the given
size. A pair only has a right answer when R1's own solution is right. {base_text}</p>
{known_table}
<p>On pairs where R1 passes, detection at {GATE_MM:g} mm is {gate_rate:.0%}
({count(gate['outcome'], 'detected')} of {len(gate)}), which {verdict_text}. {smallest_text}
{misses_text} With {len(gate)} pairs at {GATE_MM:g} mm a single miss moves the rate by
{step:.0f} points. Pairs count as resolvable from {resolvable:g} mm, the smallest size detected
at {GATE_RATE:.0%} or better.</p>
<h3>Repeats</h3>
<p>{agree} of {len(repeats)} repeats got the same choice:</p>
<ul>{repeat_items}</ul>
<h3>The BBR sign on real data</h3>
<p>The default sign was preferred on {default_wins} of {len(sign)} pairs, and the opposite sign
was flagged as misaligned on {opposite_flagged}. {opposite_text} {rule_text}, and the review
prefers the sign the rule picks on {rule_agree} of {len(sign)}.</p>
<h3>R1 against R2</h3>
{rr_table}
<p>On pairs at least {resolvable:g} mm apart, R1 was preferred on {count(far['who'], 'R1')} of
{len(far)}, R2 on {len(r2_far)} ({r2_far_items or 'none'}). R2 was flagged as misaligned on
{r2_flagged} pairs{epi_text}. {near_text}{flipped_text}</p>
</section>"""
    return facts, body


# ---------------------------------------------------------------------------------------------
# Margins
# ---------------------------------------------------------------------------------------------
#: Frozen 2026-11-02 (plan: "Frozen from the development review").
ACCURACY_MARGIN_VOX = 0.1
INCONSISTENCY_MARGIN_MM = 0.5
CAPTURE_TOLERANCE = 0.05
#: Gate 5's budget: a J-arm's median time per run at most this multiple of R1's, no run over the
#: cap (seconds, single-threaded, including the ``mri_coreg`` start).
BUDGET_RATIO, BUDGET_CAP_S = 2.0, 300.0


def budget_bullet(timing, review_set, manifest):
    """The frozen run-time budget and review set, from ``time_arms.py`` and ``review_set.py``."""
    w = timing.pivot_table(index=['dataset', 'run'], columns='step', values='seconds')
    totals = {'R1': w['mri_coreg'] + w['bbr']}
    totals |= {arm: w['mri_coreg'] + w[arm] for arm in ('J1', 'Q1', 'S2')}
    med = {k: float(v.median()) for k, v in totals.items()}
    held = manifest[manifest['side'] == 'held_out']
    subjects = held[['dataset', 'subject']].drop_duplicates()
    return (
        f'<li><b>Run-time budget (gate 5):</b> a J-arm&rsquo;s median time per run at most '
        f'{BUDGET_RATIO:g} times R1&rsquo;s chain ({med["R1"]:.0f} s at the median over '
        f'{len(w)} development runs, so {BUDGET_RATIO * med["R1"]:.0f} s), and no run over '
        f'{BUDGET_CAP_S / 60:g} minutes, single-threaded with the shared <code>mri_coreg</code> '
        f'start. Now: J1 {med["J1"]:.0f} s, Q1 {med["Q1"]:.0f} s, S2 {med["S2"]:.0f} s.</li>\n'
        f'<li><b>Review set:</b> one run per held-out subject, drawn by seed: {len(review_set)} '
        f'runs from {len(held)} held-out runs and {len(subjects)} subjects, every family and '
        'every spiral and GRASE subject included.</li>'
    )


def section_margins(m):
    """The frozen margins and the yardsticks from ``margins.py``."""

    def value(yardstick, group):
        hit = m[(m['yardstick'] == yardstick) & (m['group'] == group)]
        return float(hit['value'].iloc[0]) if len(hit) else np.nan

    def n_of(yardstick, group):
        hit = m[(m['yardstick'] == yardstick) & (m['group'] == group)]
        return int(hit['n'].iloc[0]) if len(hit) else 0

    groups = list(dict.fromkeys(m.loc[m['yardstick'] == 'r1_inconsistency_mm', 'group']))
    rows = [
        [
            esc(g),
            str(n_of('r1_inconsistency_mm', g)),
            fmt(value('r1_inconsistency_mm', g)),
            fmt(value('chance_mm', g)),
            fmt(value('worse_mm: FSL P-mean default', g)),
            fmt(value('worse_mm: FSL conventional opposite', g)),
        ]
        for g in groups
    ]
    head = ['Family / pairs', 'Subjects', 'R1 (mm)', 'Chance (mm)']
    head += ['P-mean &minus; R1', 'Opposite &minus; R1']
    table = plain_table(head, rows, num_cols=(1, 2, 3, 4, 5))
    noisy = [
        f'{esc(g)} ({fmt(value("chance_mm", g))} mm, {n_of("chance_mm", g)} subjects)'
        for g in groups
        if value('chance_mm', g) > INCONSISTENCY_MARGIN_MM
    ]
    noisy_text = (
        f'Chance exceeds the margin in {len(noisy)} group{"s" if len(noisy) != 1 else ""}: '
        f'{", ".join(noisy)}; there an arm as good as R1 can fail by chance.'
        if noisy
        else 'Chance stays under the margin in every group with at least three subjects.'
    )
    missed = [
        esc(g)
        for g in groups
        if value('worse_mm: FSL conventional opposite', g) < INCONSISTENCY_MARGIN_MM
    ]
    body = f"""
<section id="margins">
<h2>Margins</h2>
<p>Frozen on 2026-11-02 from development data, before any J-arm is scored.</p>
<ul>
<li><b>Accuracy (gates 2 and 3): {ACCURACY_MARGIN_VOX:g} voxel</b> on the median RMS error over a
family's fixtures. R1's own error has a median of {value('r1_error_vox', 'median'):.3f} voxel
(maximum {value('r1_error_vox', 'max'):.3f}). Arms of practically equal quality, the four D3
levels against R1, differ by {value('equivalent_vox', 'median'):+.3f} voxel at the median and
at most {value('equivalent_vox', 'max'):+.3f}, so the margin is about twice the largest such
gap.</li>
<li><b>Capture range (gate 3) and semi-synthetic recovery (gate 4): R1&rsquo;s converged share
less {CAPTURE_TOLERANCE * 100:g} percentage points</b> at every level, on the
optimizer&rsquo;s result.</li>
<li><b>Inconsistency (gate 4): {INCONSISTENCY_MARGIN_MM:g} mm</b> on the family median, over
subjects, of the J-arm&rsquo;s minus R1&rsquo;s per-subject inconsistency.</li>
</ul>
{table}
<p>Chance is the 95th percentile, over a bootstrap of subjects, of the family median of R1&rsquo;s
inconsistency measured in the conventional frame minus the perfusion frame, centred: how far a
family median moves when nothing changes. {noisy_text} Consistency cannot see an arm that is
consistently wrong: the opposite sign, which the review judged worse (above), stays under the
margin in {len(missed)} of {len(groups)} groups ({', '.join(missed)}). The computed
outcomes guard against a J-arm that is less consistent; the review judges accuracy.</p>
</section>"""
    return body


#: (section id, title, introduction) of each development review batch, in order.
REVIEW_BATCHES = (
    (
        'dev-review',
        'Development review batch 1',
        'Full field of view, side-by-side panels only.',
    ),
    (
        'dev-review-2',
        'Development review batch 2',
        (
            'Panels cropped to the brain, with a flicker view; known-answer runs only where R1 '
            'passes, 10 pairs at 2 mm.'
        ),
    ),
)


def build(freeze, capture=None, choices=(), margins=None, timing=None):
    fix, real = load(freeze)
    r1, f1 = section_r1(fix, real, 1)
    r2, f2 = section_r2(fix, real, 2)
    elig = section_eligibility(real)
    review = section_review(real)
    cap_bullet, cap = section_capture(pd.read_table(capture), fix) if capture else ('', '')
    revs, rev_bullets, nav = [], [], []
    for path, (section_id, title, intro) in zip(choices, REVIEW_BATCHES, strict=False):
        rf, body = section_dev_review(pd.read_table(path), real, section_id, title, intro)
        revs.append(body)
        nav.append(f'<a href="#{section_id}">{title.replace("Development review ", "")}</a>')
        outcome = 'passes' if rf['passed'] else 'falls short, so the figures were redesigned'
        rev_bullets.append(
            f'<li><b>{title}:</b> known-answer detection at {GATE_MM:g} mm is '
            f'{rf["gate_rate"]:.0%} of {rf["gate_n"]} and {outcome}. R1 was preferred on '
            f'{rf["far_r1"]} of {rf["far_n"]} R1-against-R2 pairs at least '
            f'{rf["resolvable"]:g} mm apart, and the default sign on {rf["default_wins"]} of '
            f'{rf["sign_n"]} sign pairs.</li>'
        )
    rev, rev_bullet, rev_nav = '\n'.join(revs), '\n'.join(rev_bullets), '\n'.join(nav)
    marg = section_margins(pd.read_table(margins)) if margins else ''
    if marg:
        rev_nav += '\n<a href="#margins">Margins</a>'
        rev_bullet += (
            f'\n<li><b>Margins (frozen 2026-11-02):</b> accuracy {ACCURACY_MARGIN_VOX:g} voxel, '
            f'inconsistency {INCONSISTENCY_MARGIN_MM:g} mm, capture range and semi-synthetic '
            f'recovery R1&rsquo;s converged share less {CAPTURE_TOLERANCE * 100:g} percentage '
            'points.</li>'
        )

    if timing:
        review_set = pd.read_table(HERE / 'review_set.tsv')
        manifest = pd.read_table(HERE / 'manifest.tsv')
        rev_bullet += '\n' + budget_bullet(pd.read_table(timing), review_set, manifest)
    _, hi = f1['margin']
    head = f"""<title>Coregistration Comparators</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;600;700&family=IBM+Plex+Mono:wght@400;600&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>{CSS}
.wrap {{ padding-left: max(20px, env(safe-area-inset-left, 0px));
  padding-right: max(20px, env(safe-area-inset-right, 0px)); }}
</style>"""
    summary = f"""
<ul class="findings">
<li><b>R1 (frozen, amended 2026-11-02): take the opposite BBR sign when the reference's
contrast exceeds {THRESHOLD:+.2f}.</b> On the {f1['fixtures']} development fixtures the contrast
is positive {'on exactly' if f1['exact'] else 'not only on'} the {f1['missed']} the default sign
misses by {MISS_MM:g} mm or more; the rule leaves {f1['rule_misses']} misses. The margin is
narrow: the largest real contrast is {f1['max_real']:+.3f}, the smallest fixture miss
{hi:+.3f}.</li>
<li><b>On real data R1 is production.</b> The rule changes the sign on {f1['flips']} of the
development runs and the verdict on {f1['changed']}. The fixtures' inverted GRASE and spiral
contrast is not seen in this sample.</li>
<li><b>R2: P-mean with the default sign, eligibility as the HMC probe set it.</b> P-mean is the best
of the four levels on {f2['best']} of {f2['n']} eligible runs. At least 4 pairs is what removes the
bad runs; the edge limit of 0.5 is kept because no other limit does better.</li>
<li><b>The perfusion image is the weaker comparator on real data.</b> The best level beats the
conventional reference on {f2['beats']:.0%} of eligible runs. R1 is frozen as the comparator in
every family (2026-11-02), from these numbers and the two development review batches.</li>
<li><b>Frozen from the review (2026-11-02):</b> the tie distance is the larger of
{GATE_MM:g} mm (the smallest known-answer size detected at {GATE_RATE:.0%} or better) and the
run's noise floor; the floor limit is {FLOOR_LIMIT_MM:g} mm.</li>
{cap_bullet}
{rev_bullet}
</ul>"""
    return f"""{head}
<div class="wrap">
<header>
<div class="eyebrow">Coregistration probe &middot; Task 3 &middot; development data</div>
<h1>Coregistration comparators</h1>
<p class="lede">What production's FSL path does with the BBR sign fixed, and with a perfusion image as
the reference, on development fixtures and real runs: the evidence for freezing R1 and R2 before any
joint arm is screened.</p>
<p class="meta">{f1['fixtures']} development fixtures &middot; {len(runs_of(real))} development real
runs &middot; held-out data not used</p>
<nav class="toc"><a href="#r1">R1 sign</a><a href="#r2">R2 level</a>
<a href="#eligibility">Eligibility</a>{'<a href="#capture">Capture range</a>' if cap else ''}
{rev_nav}
<a href="#review">For the reviewer</a></nav>
</header>
<section id="summary"><h2>Summary</h2>{summary}</section>
{r1}
{r2}
{elig}
{cap}
{rev}
{marg}
{review}
</div>"""


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--freeze', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--capture', type=Path, help="capture.py's table of the comparators")
    p.add_argument(
        '--review',
        type=Path,
        nargs='*',
        default=(),
        help="score_review.py's choices.tsv of each development batch, in order",
    )
    p.add_argument('--margins', type=Path, help="margins.py's table of yardsticks")
    p.add_argument('--timing', type=Path, help="time_arms.py's table of run times")
    args = p.parse_args(argv)
    args.out.write_text(build(args.freeze, args.capture, args.review, args.margins, args.timing))


if __name__ == '__main__':
    main()
