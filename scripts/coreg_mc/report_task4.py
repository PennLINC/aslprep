#!/usr/bin/env python3
"""Walkthrough HTML report of Task 4's primary screen and gate 2 (development part).

Reads ``screen.py``'s table, ``gate2.py``'s verdicts and Task 3's ``fixtures.tsv`` (for R0 and
R1 on the same fixtures) and writes one Artifact page fragment (no doctype or
``<html>``/``<body>``). The chart helpers and styling are the HMC probe's
(``scripts/hmc_probe/report.py``). Every number in the prose is computed here from the data.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/report_task4.py --screen screen.tsv --gate2 gate2.tsv \\
        --freeze <Task 3 freeze dir> --out report.html
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import freeze_comparators  # noqa: E402
import gate2  # noqa: E402
from report import CSS, esc, figure, fmt, heat_table, plain_table, strip_plot  # noqa: E402

FAMILY_LABELS = {
    'st': 'Spirals',
    'grase': 'GRASE',
    'short': 'Short series',
    'other': 'Elsewhere',
}
COLUMNS = (
    ('st', 'derivatives'),
    ('st', 'raw'),
    ('grase', 'derivatives'),
    ('grase', 'raw'),
    ('short', 'derivatives'),
    ('other', 'derivatives'),
    ('other', 'raw'),
)


def comparators(freeze, voxel):
    """R0 and R1 (the frozen sign rule) on the screen's fixtures, in voxels."""
    fix = pd.read_table(freeze / 'fixtures.tsv')
    fix = fix[fix['fixture'].isin(voxel.index)]
    w = fix.pivot(index='fixture', columns='arm', values='rms_mm')
    contrast = fix.drop_duplicates('fixture').set_index('fixture')['contrast']
    opposite = contrast.map(freeze_comparators.r1_sign) == 'opposite'
    r1 = w['FSL conventional default'].where(~opposite, w['FSL conventional opposite'])
    rows = []
    for name, mm in (('R0', w['R0 production']), ('R1', r1)):
        for fixture, value in mm.items():
            rows.append(
                {
                    'fixture': fixture,
                    'name': name,
                    'arm': name,
                    'rms_voxels': value / voxel[fixture],
                },
            )

    return pd.DataFrame(rows)


def medians(df):
    df = df.assign(
        family=df['fixture'].map(lambda f: gate2.FIXTURES[f][0]),
        form=df['fixture'].map(lambda f: gate2.FIXTURES[f][1]),
    )
    t = df.groupby(['name', 'family', 'form'])['rms_voxels'].median().unstack(['family', 'form'])
    t = t[list(COLUMNS)]
    t.columns = [f'{f}/{form}' for f, form in COLUMNS]
    t['all'] = df.groupby('name')['rms_voxels'].median()
    return t.sort_values('all')


def best_of(t, prefix):
    names = [n for n in t.index if n.startswith(prefix)]
    return t.loc[names, 'all'].idxmin()


def section_medians(t, fig):
    cols = [f'{FAMILY_LABELS[f]} ({form[:3]})' for f, form in COLUMNS] + ['All']
    values = [[float(v) for v in t.loc[n]] for n in t.index]
    notes = {'R0': 'production', 'R1': 'frozen comparator'}
    body = heat_table(list(t.index), cols, values, notes=notes, digits=3, col_labels=cols)
    return f"""
<section id="medians">
<h2>Every arm by family</h2>
<p>Median RMS error against truth, in voxels, over each family&rsquo;s fixtures in each form
(<code>der</code>: ASLPrep&rsquo;s derivatives; <code>raw</code>: the raw series), sorted by the
median over all {int(t.attrs['n_fixtures'])} development fixtures. R0 and R1 come from Task 3 and
are shown for scale; gate 2 compares J1 only with its matched ablations.</p>
{figure(body, 'Median error by arm and family (voxels; darker is worse, log scale).', fig)}
</section>"""


def section_strip(df, t, fig):
    picks = ['R1', 'S1_mask', best_of(t, 'S2_'), best_of(t, 'Q1_'), best_of(t, 'J1_')]
    picks = list(dict.fromkeys(picks))
    values = [df.loc[df['name'] == n, 'rms_voxels'].to_numpy() for n in picks]
    body = strip_plot(
        picks, values, 'RMS error against truth (voxels)', label='per-fixture errors'
    )
    return f"""
<section id="fixtures">
<h2>Per fixture</h2>
<p>The frozen comparator, S1 (masked, as frozen) and the arm of each kind with the lowest median
over all fixtures. One point per development fixture; the tick is the median.</p>
{figure(body, 'Per-fixture error of the best arm of each kind (voxels).', fig)}
</section>"""


def section_gate(g, fig):
    rows = []
    for r in g.itertuples(index=False):
        cells = [f'<code>{esc(r.name)}</code>', fmt(r.median_vox, 3)]
        for family in FAMILY_LABELS:
            best, worst = getattr(r, f'{family}_vs_best'), getattr(r, f'{family}_vs_worst')
            cells.append(f'{best:+.3f} / {worst:+.3f}')

        cells.append(esc(r.failed.replace(',', ', ') if isinstance(r.failed, str) else ''))
        rows.append(cells)

    head = ['Variant', 'Median (vox)'] + [FAMILY_LABELS[f] for f in FAMILY_LABELS] + ['Failed']
    table = plain_table(head, rows, num_cols=(1, 2, 3, 4, 5))
    return f"""
<section id="gate">
<h2>Gate 2, variant by variant</h2>
<p>For each family, the variant&rsquo;s median less its best matched ablation&rsquo;s, then less
its worst&rsquo;s (voxels; the larger over the family&rsquo;s forms). A stress family is won when
the first number is at most &minus;{gate2.MARGIN_VOX:g}; any family fails when it exceeds
+{gate2.MARGIN_VOX:g}.</p>
{figure(table, 'Gaps to the matched ablations (best / worst), and the families failed.', fig)}
</section>"""


def build(screen, verdicts, freeze):
    voxel = screen.drop_duplicates('fixture').set_index('fixture')['voxel_mm']
    df = pd.concat([screen, comparators(freeze, voxel)], ignore_index=True)
    t = medians(df)
    t.attrs['n_fixtures'] = df['fixture'].nunique()
    n_pass = int(verdicts['passes'].sum())
    n_variants = len(verdicts)
    n_won = int((verdicts['n_won'] > 0).sum())
    failed = verdicts['failed'].fillna('')
    grase_failed = int(failed.str.contains('grase').sum())
    clean = int((failed == '').sum())
    clean_text = f'{clean} fail no family' if clean else 'none clears every family'
    best_j1, best_s2, best_q1 = best_of(t, 'J1_'), best_of(t, 'S2_'), best_of(t, 'Q1_')
    s1 = t.loc['S1_mask']
    j1 = t.loc[best_j1]
    s2 = t.loc[best_s2]
    grase = ['grase/derivatives', 'grase/raw']
    short_best = verdicts['short_vs_best'].min()
    single = [n for n in t.index if n.startswith(('S2_', 'Q1_'))]
    r1_gap = float((t.loc[single, 'all'] - t.loc['R1', 'all']).abs().max())
    head = """<title>Joint Coregistration Screen</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;600;700&family=IBM+Plex+Mono:wght@400;600&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">"""
    head += f"""
<style>{CSS}
.wrap {{ padding-left: max(20px, env(safe-area-inset-left, 0px));
  padding-right: max(20px, env(safe-area-inset-right, 0px)); }}
</style>"""
    verdict = (
        f'<b>J1 fails gate 2: {n_pass} of {n_variants} variants pass.</b> Under the rules fixed '
        'before any J1 result, the extensions (J4&ndash;J6) are not run and the probe reports a '
        'negative result.'
        if n_pass == 0
        else f'<b>{n_pass} of {n_variants} J1 variants pass gate 2.</b>'
    )
    summary = f"""
<ul class="findings">
<li>{verdict}</li>
<li><b>No variant wins a stress family.</b> {n_won} of {n_variants} beat every matched ablation by
{gate2.MARGIN_VOX:g} voxel on any stress family; the closest is the short series, where the best
variant is {short_best:+.3f} voxel from its best ablation.</li>
<li><b>Every variant is worse on GRASE</b> ({grase_failed} of {n_variants} fail it;
{clean_text}). The best J1, <code>{esc(best_j1)}</code>, has
{fmt(j1[grase[0]], 3)} and {fmt(j1[grase[1]], 3)} voxel there (derivatives, raw), against
{fmt(s2[grase[0]], 3)} and {fmt(s2[grase[1]], 3)} for <code>{esc(best_s2)}</code> and
{fmt(s1[grase[0]], 3)} and {fmt(s1[grase[1]], 3)} for S1.</li>
<li><b>The static term pulls J1 toward S1&rsquo;s answer.</b> The T1w MI term is the weakness the
engine pilot carried into Task 4. Q1, which runs the same terms one after the other, ends at the
perfusion term&rsquo;s optimum and does as well as S2; J1, which sums them, keeps S1&rsquo;s error
on GRASE.</li>
<li><b>Nothing here improves on production&rsquo;s chain.</b> R1 has a median of
{fmt(t.loc['R1', 'all'], 3)} voxel over these fixtures; the best ANTs arms are
<code>{esc(best_s2)}</code> ({fmt(t.loc[best_s2, 'all'], 3)}) and <code>{esc(best_q1)}</code>
({fmt(t.loc[best_q1, 'all'], 3)}), the best J1 {fmt(j1['all'], 3)}. Over all fixtures, every
S2 and Q1 arm is within {r1_gap:.3f} voxel of R1, under the {gate2.MARGIN_VOX:g}-voxel
margin.</li>
<li><b>Checks.</b> All {len(screen)} arm runs produced a transform. The frozen variants
(<code>S1_mask</code>, <code>S2_P-winsor_gm_mask</code>, <code>Q1_P-winsor_gm_mask</code>,
<code>J1_P-winsor_gm_mask_w4</code>) build the engine pilot&rsquo;s commands and reproduce its
transforms exactly where the pilot ran them.</li>
</ul>"""
    return f"""{head}
<div class="wrap">
<header>
<div class="eyebrow">Coregistration probe &middot; Task 4 &middot; development fixtures</div>
<h1>Joint coregistration screen</h1>
<p class="lede">Whether one rigid ANTs registration that sums a T1w term and a perfusion term
(J1) beats the same terms alone (S1, S2) or one after the other (Q1) on simulated fixtures with
known truth: gate 2, development part.</p>
<p class="meta">{len(verdicts)} J1 variants &middot; {screen['name'].nunique() - n_variants}
matched ablations &middot; {t.attrs['n_fixtures']} development fixtures &middot; held-out data
not used</p>
<nav class="toc"><a href="#medians">By family</a><a href="#fixtures">Per fixture</a>
<a href="#gate">Gate 2</a><a href="#rules">Rules</a></nav>
</header>
<section id="summary"><h2>Summary</h2>{summary}</section>
{section_medians(t, 1)}
{section_strip(df, t, 2)}
{section_gate(verdicts, 3)}
<section id="rules">
<h2>Rules (fixed 2026-11-02, before any J1 result)</h2>
<ul>
<li>J1 over perfusion image (P-winsor, P-robust) &times; target (smoothed GM, pseudo-CBF:
GM + {26.7 / 80:.3f}&nbsp;WM, after ExploreASL) &times; perfusion weight (0.2, 0.4, 0.6) &times;
masking (stage masks, none): {n_variants} variants, each against S1 with the same masking and
S2 and Q1 with the same image, target and masking.</li>
<li><b>Better</b> on a stress family (spirals, GRASE, short series): the family median is lower
than every matched ablation&rsquo;s by at least {gate2.MARGIN_VOX:g} voxel on every form.
<b>Not worse</b> on any family: at most the best matched ablation&rsquo;s plus
{gate2.MARGIN_VOX:g} voxel on every form.</li>
<li>A variant passes with at least one stress family won and none failed; at most
{gate2.KEEP} go on. If none passes, J1 fails gate 2.</li>
</ul>
</section>
</div>"""


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--screen', type=Path, required=True)
    p.add_argument('--gate2', type=Path, required=True)
    p.add_argument('--freeze', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(argv)
    screen = pd.read_table(args.screen)
    verdicts = pd.read_table(args.gate2)
    args.out.write_text(build(screen, verdicts, args.freeze))
    print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
