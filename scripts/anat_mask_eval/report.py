#!/usr/bin/env python3
"""Build the walkthrough report (one HTML file) from ``evaluate.py``'s results.

Every number in the text is computed here from the result tables. The only hand-entered
content is ``VERDICTS``: what the montage of each strongly disagreeing reference shows.

Usage::

    python report.py --out-dir <anat-mask-eval dir> --html report.html
"""

import argparse
import base64
import io
import json
import re
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd
from svgplot import METHOD_STYLE, bands, strip

#: Below this Dice between the anatomical and SynthStrip masks, a reference is a disagreement.
DISAGREE = 0.90

#: Montage reading of every disagreement: which mask is wrong, and how.
VERDICTS = {
    'ds004856__sub-3131_ses-wave1_acq-ASL_run-1': (
        'coreg',
        'Anatomical mask displaced: coregistration failed (BBR fell back to a wrong initialisation)',
    ),
    'ds004856__sub-3131_ses-wave2_acq-ASL_run-1': (
        'coreg',
        'Anatomical mask displaced: coregistration failed (BBR fell back)',
    ),
    'ds004856__sub-3131_ses-wave3_acq-ASL_run-1': (
        'coreg',
        'Anatomical mask displaced: coregistration failed (BBR fell back)',
    ),
    'ds008719__sub-050_echo-1': (
        'coreg',
        'Anatomical mask shrunk and shifted inferiorly: coregistration wrong, with no fallback',
    ),
    'lsa-aslprep-debug__sub-761_ses-1': (
        'synthstrip',
        'SynthStrip mask cut by a flat plane through the brain',
    ),
    'tarricone__sub-Sub004': ('synthstrip', 'SynthStrip mask with large holes and missing cortex'),
    'morys_bad_brainmask__sub-MNI150003': (
        'synthstrip',
        'SynthStrip mask cut by a flat plane through the brain',
    ),
    'jsr-aslprep-debug__sub-lsa116_ses-03_acq-3dPCASL': (
        'synthstrip',
        'SynthStrip mask cut by a flat plane through the posterior brain',
    ),
}
RIGHT = {'coreg': 'synthstrip', 'synthstrip': 'anat_pv50'}


def _img(path, alt):
    """A montage embedded as JPEG (a third of the PNG's size, for a page under 16 MB)."""
    from PIL import Image

    buf = io.BytesIO()
    Image.open(path).convert('RGB').save(buf, 'JPEG', quality=85)
    data = base64.b64encode(buf.getvalue()).decode()
    return f'<img src="data:image/jpeg;base64,{data}" alt="{escape(alt)}" loading="lazy">'


def _f(v, nd=3):
    return '—' if v is None or not np.isfinite(v) else f'{v:.{nd}f}'


def short(name):
    return name.replace('__', ' · ')


def load(out):
    r = out / 'results'
    masks = pd.read_csv(r / 'real_masks.tsv', sep='\t')
    runs = pd.read_csv(r / 'real_runs.tsv', sep='\t')
    perturb = pd.read_csv(r / 'real_perturb.tsv', sep='\t')
    sim = pd.read_csv(r / 'sim_masks.tsv', sep='\t')
    problems = pd.read_csv(r / 'problems.tsv', sep='\t')
    # One row per distinct reference image: runs sharing a reference share every mask.
    first = runs.drop_duplicates('ref_sha1').set_index('name')
    masks = masks[masks['name'].isin(first.index)]
    perturb = perturb[perturb['name'].isin(first.index)]
    return masks, first, perturb, sim, problems, runs


def table(header, rows):
    head = ''.join(f'<th>{escape(h)}</th>' for h in header)
    body = ''.join('<tr>' + ''.join(rows_i) + '</tr>' for rows_i in rows)
    return f'<div class="tablewrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def build(out, html_path):
    masks, runs, perturb, sim, problems, all_runs = load(out)

    def wide(col):
        return masks.pivot(index='name', columns='method', values=col).loc[runs.index]

    d_ss, vol = wide('dice_synthstrip'), wide('vol_ml')
    d_anat = wide('dice_anat_pv50')
    disagree = d_ss.index[d_ss['anat_pv50'] < DISAGREE]
    ok = d_ss.index.difference(disagree)
    missing = [n for n in disagree if n not in VERDICTS]
    if missing:
        raise SystemExit(f'Look at these montages and add verdicts: {missing}')

    def med(s):
        return float(np.nanmedian(s))

    # ---- Ground truth ---------------------------------------------------------------------
    sim_w = sim.pivot(index='name', columns='method', values='dice_truth')
    sim_methods = [
        'synthstrip',
        'legacy',
        'anat_nn',
        'anat_pv50',
        'anat_pv10',
        'hybrid',
        'hybrid_vol',
        'anat_truexfm',
    ]
    sim_err = sim.groupby('name')['coreg_err_mm'].first()
    bad_sim = sim_err.index[sim_err > 3]
    good_sim = sim_err.index.difference(bad_sim)
    chart_sim = strip(
        [
            (
                METHOD_STYLE[m][0],
                METHOD_STYLE[m][1],
                [(short(n), sim_w.loc[n, m]) for n in sim_w.index],
            )
            for m in sim_methods
        ],
        'Dice with the true brain (one dot per simulated run)',
        0.85,
        1.0,
    )

    def sim_rows():
        rows = []
        for m in sim_methods:
            sub = sim[(sim['method'] == m) & sim['name'].isin(good_sim)]
            rows.append(
                [
                    f'<th scope="row">{METHOD_STYLE[m][0]}</th>',
                    f'<td>{_f(sub["dice_truth"].median())}</td>',
                    f'<td>{_f(sub["false_pos_ml"].median(), 0)}</td>',
                    f'<td>{_f(sub["false_neg_ml"].median(), 0)}</td>',
                    f'<td>{_f(sub["gm_kept"].median())}</td>',
                ]
                + [f'<td>{_f(sim_w.loc[n, m])}</td>' for n in bad_sim]
            )
        return rows

    sim_table = table(
        [
            'Mask',
            'Dice with truth',
            'Non-brain included (mL)',
            'Brain missed (mL)',
            'Grey matter kept',
        ]
        + [f'{short(n)} (coreg off {sim_err[n]:.0f} mm): Dice' for n in bad_sim],
        sim_rows(),
    )

    # ---- Real-data agreement --------------------------------------------------------------
    chart_dice = strip(
        [
            (
                f'{METHOD_STYLE[m][0]} vs SynthStrip',
                METHOD_STYLE[m][1],
                [(short(n), d_ss.loc[n, m]) for n in d_ss.index],
            )
            for m in ('legacy', 'anat_nn', 'anat_pv50')
        ],
        'Dice with the SynthStrip mask (one dot per reference image)',
        0.4,
        1.0,
        ref=DISAGREE,
    )
    ratio = vol.div(vol['synthstrip'], axis=0)
    chart_vol = strip(
        [
            (METHOD_STYLE[m][0], METHOD_STYLE[m][1], [(short(n), ratio.loc[n, m]) for n in ok])
            for m in ('legacy', 'anat_pv50')
        ],
        'Mask volume ÷ SynthStrip mask volume (agreeing references)',
        0.7,
        1.5,
        fmt='{:.2f}',
        ref=1.0,
    )

    fam = runs['family']

    def family_rows():
        g = pd.DataFrame(
            {
                'family': fam,
                'anat': d_ss['anat_pv50'],
                'legacy': d_ss['legacy'],
                'sdc': runs['sdc'],
            }
        )
        rows = []
        for f, sub in g.groupby('family'):
            rows.append(
                [
                    f'<th scope="row">{escape(f)}</th>',
                    f'<td>{len(sub)}</td>',
                    f'<td>{int(sub["sdc"].sum())}</td>',
                    f'<td>{_f(sub["anat"].median())}</td>',
                    f'<td>{_f(sub["anat"].min())}</td>',
                    f'<td>{_f(sub["legacy"].median())}</td>',
                ]
            )
        return rows

    family_table = table(
        [
            'Acquisition family',
            'References',
            'Distortion-corrected',
            'Anatomical vs SynthStrip: median Dice',
            'lowest',
            'Legacy vs SynthStrip: median Dice',
        ],
        family_rows(),
    )

    # ---- Registration sensitivity ---------------------------------------------------------
    shifts = sorted(perturb['shift_mm'].unique())
    by_shift = {0.0: np.ones(1)} | {
        s: perturb.loc[perturb['shift_mm'] == s, 'dice'].to_numpy() for s in shifts
    }
    chart_perturb = bands(
        [('Anatomical mask', 'var(--s3)', by_shift)],
        [0.0, *shifts],
        'Rigid shift added to the coregistration (mm)',
        'Dice with the unshifted mask',
        0.88,
        1.0,
    )
    curve_x = np.array([0.0, *shifts])
    curve_y = np.array([1.0, *[np.median(by_shift[s]) for s in shifts]])
    equiv_mm = float(np.interp(-med(d_ss.loc[ok, 'anat_pv50']), -curve_y, curve_x))

    # ---- Disagreements --------------------------------------------------------------------
    def disagree_rows():
        rows = []
        for n in sorted(disagree, key=lambda n: d_ss.loc[n, 'anat_pv50']):
            kind, text = VERDICTS[n]
            right = RIGHT[kind]
            picks = []
            for g in ('hybrid', 'hybrid_vol'):
                chose = 'anat_pv50' if d_anat.loc[n, g] == 1 else 'synthstrip'
                picks.append(
                    '<td class="ok">right</td>' if chose == right else '<td class="bad">wrong</td>'
                )
            rows.append(
                [
                    f'<th scope="row">{escape(short(n))}</th>',
                    f'<td>{escape(fam[n])}</td>',
                    f'<td>{_f(d_ss.loc[n, "anat_pv50"])}</td>',
                    f'<td>{_f(ratio.loc[n, "anat_pv50"] ** -1, 2)}</td>',
                    f'<td>{_f(runs.loc[n, "alt_dice_min"])}</td>',
                    f'<td>{"yes" if runs.loc[n, "coreg_fallback"] else "no"}</td>',
                    f'<td>{escape(text)}</td>',
                    *picks,
                ]
            )
        return rows

    disagree_table = table(
        [
            'Reference',
            'Family',
            'Dice anatomical vs SynthStrip',
            'SynthStrip ÷ anatomical volume',
            'Lowest Dice over other registrations',
            'BBR fallback',
            'What the montage shows',
            'Dice guard',
            'Dice + volume guard',
        ],
        disagree_rows(),
    )
    guard_right = {
        g: sum(
            ('anat_pv50' if d_anat.loc[n, g] == 1 else 'synthstrip') == RIGHT[VERDICTS[n][0]]
            for n in disagree
        )
        for g in ('hybrid', 'hybrid_vol')
    }
    n_coreg = sum(VERDICTS[n][0] == 'coreg' for n in disagree)
    n_ssfail = sum(VERDICTS[n][0] == 'synthstrip' for n in disagree)

    # ---- Perfusion edges ------------------------------------------------------------------
    perf_charts = {}
    for col, axis, lo, hi in (
        (
            'outer_ratio',
            'Mean perfusion in the one-voxel shell outside the mask ÷ mean inside',
            -0.4,
            1.0,
        ),
        ('inner_ratio', "Mean perfusion in the mask's own edge voxels ÷ mean inside", -0.4, 1.2),
    ):
        w = wide(col).loc[ok]
        perf_charts[col] = strip(
            [
                (
                    METHOD_STYLE[m][0],
                    METHOD_STYLE[m][1],
                    [(short(n), w.loc[n, m]) for n in w.index],
                )
                for m in ('synthstrip', 'legacy', 'anat_pv50')
            ],
            axis,
            lo,
            hi,
            fmt='{:.1f}',
        )

    def summary_rows():
        cols = [
            ('vol_ml', 0),
            ('dice_synthstrip', 3),
            ('outer_ratio', 2),
            ('inner_ratio', 2),
            ('perf_out_frac', 2),
        ]
        rows = []
        for m in ('synthstrip', 'legacy', 'anat_nn', 'anat_pv50', 'anat_pv10'):
            sub = masks[(masks['method'] == m) & masks['name'].isin(ok)]
            rows.append(
                [f'<th scope="row">{METHOD_STYLE[m][0]}</th>']
                + [f'<td>{_f(sub[c].median(), nd)}</td>' for c, nd in cols]
            )
        return rows

    summary_table = table(
        [
            'Mask (medians over agreeing references)',
            'Volume (mL)',
            'Dice vs SynthStrip',
            'Outside ÷ inside',
            'Edge ÷ inside',
            'Perfusion signal outside the mask',
        ],
        summary_rows(),
    )

    # ---- Galleries ------------------------------------------------------------------------
    fig_dir = out / 'figures'

    def figure(n, cap):
        p = fig_dir / f'{n}.png'
        return (
            f'<figure>{_img(p, short(n))}<figcaption>{escape(cap)}</figcaption></figure>'
            if p.exists()
            else ''
        )

    disagree_gallery = '\n'.join(
        figure(n, f'{short(n)}: {VERDICTS[n][1]} (Dice {_f(d_ss.loc[n, "anat_pv50"])})')
        for n in sorted(disagree, key=lambda n: d_ss.loc[n, 'anat_pv50'])
    )
    sim_gallery = '\n'.join(
        figure(
            n,
            f'{short(n)}: Dice with truth, SynthStrip {_f(sim_w.loc[n, "synthstrip"])}, '
            f'anatomical {_f(sim_w.loc[n, "anat_pv50"])}; coregistration error {sim_err[n]:.1f} mm',
        )
        for n in sim_w.index
    )
    ok_gallery = '\n'.join(
        figure(
            n,
            f'{short(n)} ({fam[n]}): anatomical vs SynthStrip Dice {_f(d_ss.loc[n, "anat_pv50"])}',
        )
        for n in sorted(ok, key=lambda n: (fam[n], n))
    )

    sdc = runs['sdc']
    stats = {
        'n_runs': len(all_runs),
        'n_refs': len(runs),
        'n_sets': int(runs['dataset'].nunique()),
        'n_families': int(fam.nunique()),
        'n_disagree': len(disagree),
        'n_ok': len(ok),
        'n_coreg': n_coreg,
        'n_ssfail': n_ssfail,
        'n_sdc': int(sdc.sum()),
        'n_sims': int(sim['name'].nunique()),
        'n_good_sim': len(good_sim),
        'guard_dice_right': int(guard_right['hybrid']),
        'guard_vol_right': int(guard_right['hybrid_vol']),
        'anat_ss_ok': med(d_ss.loc[ok, 'anat_pv50']),
        'legacy_ss_ok': med(d_ss.loc[ok, 'legacy']),
        'nn_vs_pv': med(d_anat.loc[ok, 'anat_nn']),
        'legacy_vol_ratio': med(ratio.loc[ok, 'legacy']),
        'anat_vol_ratio': med(ratio.loc[ok, 'anat_pv50']),
        'anat_ss_sdc': med(d_ss.loc[ok.intersection(sdc.index[sdc]), 'anat_pv50']),
        'anat_ss_nosdc': med(d_ss.loc[ok.intersection(sdc.index[~sdc]), 'anat_pv50']),
        'prod_is_legacy_sdc': int(
            (d_ss.loc[sdc.index[sdc], 'production'] - d_ss.loc[sdc.index[sdc], 'legacy'])
            .abs()
            .lt(1e-9)
            .sum()
        ),
        'prod_ss_nosdc': med(d_ss.loc[sdc.index[~sdc], 'production']),
        'equiv_mm': equiv_mm,
        'p1': med(by_shift[1.0]),
        'p2': med(by_shift[2.0]),
        'p3': med(by_shift[3.0]),
        'p6': med(by_shift[6.0]),
        'sim_err_good_max': float(sim_err[good_sim].max()),
        'sim_bad_name': short(bad_sim[0]) if len(bad_sim) else '—',
        'sim_bad_err': float(sim_err[bad_sim].max()) if len(bad_sim) else np.nan,
        'sim_bad_anat': float(sim_w.loc[bad_sim, 'anat_pv50'].min()) if len(bad_sim) else np.nan,
        'sim_bad_ss': float(sim_w.loc[bad_sim, 'synthstrip'].min()) if len(bad_sim) else np.nan,
        'perf_outer_ss': med(wide('outer_ratio').loc[ok, 'synthstrip']),
        'perf_outer_anat': med(wide('outer_ratio').loc[ok, 'anat_pv50']),
        'perf_outer_leg': med(wide('outer_ratio').loc[ok, 'legacy']),
        'perf_inner_ss': med(wide('inner_ratio').loc[ok, 'synthstrip']),
        'perf_inner_anat': med(wide('inner_ratio').loc[ok, 'anat_pv50']),
        'perf_inner_leg': med(wide('inner_ratio').loc[ok, 'legacy']),
        'in_fov_min': float(runs['anat_in_fov'].min()),
    }
    for m in sim_methods:
        sub = sim[(sim['method'] == m) & sim['name'].isin(good_sim)]
        stats[f'sim_{m}'] = float(sub['dice_truth'].median())
        stats[f'simgm_{m}'] = float(sub['gm_kept'].median())
        stats[f'simfp_{m}'] = float(sub['false_pos_ml'].median())
        stats[f'simfn_{m}'] = float(sub['false_neg_ml'].median())
    (out / 'results' / 'report_stats.json').write_text(json.dumps(stats, indent=1, default=float))

    fmt = {
        k: (
            v
            if isinstance(v, (int, str))
            else _f(v, 0 if k.startswith(('simfp', 'simfn', 'sim_bad_err')) else 3)
        )
        for k, v in stats.items()
    }
    fmt['equiv_mm'] = f'{equiv_mm:.1f}'
    for k in (
        'legacy_vol_ratio',
        'anat_vol_ratio',
        'perf_outer_ss',
        'perf_outer_anat',
        'perf_outer_leg',
        'perf_inner_ss',
        'perf_inner_anat',
        'perf_inner_leg',
        'in_fov_min',
    ):
        fmt[k] = _f(stats[k], 2)
    for kind in ('coreg', 'synthstrip'):
        r = [1 / ratio.loc[n, 'anat_pv50'] for n in disagree if VERDICTS[n][0] == kind]
        fmt[f'ratio_{kind}'] = f'{min(r):.2f}–{max(r):.2f}'  # noqa: RUF001
    fmt['sim_err_good_max'] = f'{stats["sim_err_good_max"]:.1f}'
    fmt['simgm_loss_ss'] = f'{100 * (1 - stats["simgm_synthstrip"]):.0f}'
    fmt['simgm_loss_anat'] = f'{100 * (1 - stats["simgm_anat_pv50"]):.0f}'

    tpl = (Path(__file__).parent / 'report_template.html').read_text(encoding='utf-8')
    values = dict(
        chart_dice=chart_dice,
        chart_vol=chart_vol,
        chart_outer=perf_charts['outer_ratio'],
        chart_inner=perf_charts['inner_ratio'],
        chart_perturb=chart_perturb,
        chart_sim=chart_sim,
        summary_table=summary_table,
        sim_table=sim_table,
        disagree_table=disagree_table,
        family_table=family_table,
        disagree_gallery=disagree_gallery,
        ok_gallery=ok_gallery,
        sim_gallery=sim_gallery,
        problems=len(problems),
        **fmt,
    )
    # ``@@name@@`` placeholders, so the template's CSS braces need no escaping.
    html = re.sub(r'@@(\w+)@@', lambda m: str(values[m.group(1)]), tpl)
    html_path.write_text(html, encoding='utf-8')
    return stats


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--html', type=Path, required=True)
    args = parser.parse_args(argv)
    stats = build(args.out_dir, args.html)
    print(json.dumps(stats, indent=1, default=float))


if __name__ == '__main__':
    main()
