#!/usr/bin/env python3
"""Walkthrough HTML report of Task 3's comparator evidence: R1's sign rule and R2's D3 level.

Reads ``freeze_comparators.py``'s ``fixtures.tsv`` and ``real.tsv`` and writes one Artifact page
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
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

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


def rule_sign(contrast):
    """R1's proposed sign: the opposite sign when white matter is brighter than grey."""
    return 'opposite' if contrast > 0 else 'default'


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
positive contrast selects the opposite sign.</p>
<p>On the {len(f)} development fixtures, the {len(inverted)} with positive contrast
({', '.join(f'<code>{esc(n)}</code>' for n in inverted.index)},
{inverted['contrast'].min():+.3f} to {inverted['contrast'].max():+.3f})
{'are exactly' if exact else 'do not match'} the {len(missed)} the default sign misses by
{MISS_MM:g} mm or more. The other {len(upright)} sit between {upright['contrast'].min():+.3f}
and {upright['contrast'].max():+.3f}. The rule takes the default sign's {len(missed)} misses to
{rule_misses}.</p>
{chart}
{fix_table}
<p>Real data does not show the inversion. Across {len(r)} development runs the contrast is
negative in every family; only {len(flips)} runs are positive, the largest at
{r['contrast'].max():+.3f}, and on {len(flips) - len(changed)} of them both signs get the same
verdict. The default sign has {total['default']} gross failures and the opposite sign
{total['opposite']}; the rule has {total['rule']}. {shared} of the failures are runs every arm
fails (<a href="#review">below</a>), so they say nothing about the sign.</p>
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


def build(freeze):
    fix, real = load(freeze)
    r1, f1 = section_r1(fix, real, 1)
    r2, f2 = section_r2(fix, real, 2)
    elig = section_eligibility(real)
    review = section_review(real)
    lo, hi = f1['margin']
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
<li><b>R1: take the opposite BBR sign when the reference's white matter is brighter than its grey
matter.</b> On the {f1['fixtures']} development fixtures the contrast is positive
{'on exactly' if f1['exact'] else 'not only on'} the {f1['missed']} the default sign misses by
{MISS_MM:g} mm or more; the rule leaves {f1['rule_misses']} misses. The margin is narrow
({lo:+.3f} against {hi:+.3f}).</li>
<li><b>On real data R1 is production.</b> The contrast is negative in every family; the rule
changes the sign on {f1['flips']} of the development runs and the verdict on {f1['changed']}.
The fixtures' inverted GRASE and spiral contrast is not seen in this sample.</li>
<li><b>R2: P-mean with the default sign, eligibility as the HMC probe set it.</b> P-mean is the best
of the four levels on {f2['best']} of {f2['n']} eligible runs. At least 4 pairs is what removes the
bad runs; the edge limit of 0.5 is kept because no other limit does better.</li>
<li><b>The perfusion image is the weaker comparator on real data.</b> The best level beats the
conventional reference on {f2['beats']:.0%} of eligible runs. The development review batch decides
the comparator per family; these numbers point to R1 everywhere.</li>
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
<a href="#eligibility">Eligibility</a><a href="#review">For the reviewer</a></nav>
</header>
<section id="summary"><h2>Proposals</h2>{summary}</section>
{r1}
{r2}
{elig}
{review}
</div>"""


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--freeze', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(argv)
    args.out.write_text(build(args.freeze))


if __name__ == '__main__':
    main()
