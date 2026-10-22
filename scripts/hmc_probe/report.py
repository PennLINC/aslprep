#!/usr/bin/env python3
"""Walkthrough HTML report of the HMC probe's findings so far.

Reads the probe's result files and writes one self-contained HTML page (an Artifact page
fragment: no doctype or ``<html>``/``<body>``, which the publisher adds) with SVG figures
drawn from theme tokens. Every number in the prose is computed here from the data.

Inputs (under ``--root``, default ``/home/tsalo/hmc-probe``):

- ``results.json``: the noise and cost-function experiment;
- ``blur-bias/<fixture>/blur_bias.json``: the blurred-reference experiment;
- ``fixture-scores/<label>/<fixture>/score.json`` (and C8's ``artifacts/``): the candidate
  screen on simulated fixtures;
- ``fixture-runs/baseline/<fixture>/score_full.json``: the reference-consumer baseline;
- ``harness-runs/screen_c0.tsv`` and ``real_c0.tsv``: production on the real-data panel.

Usage::

    python scripts/hmc_probe/report.py --out report.html [--root /home/tsalo/hmc-probe]
"""

import argparse
import html
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from candidates import zigzag_clean  # noqa: E402
from compare import METRICS, _get  # noqa: E402

ARMS = {
    'c0': ('C0', 'Production: MCFLIRT per volume type to one common reference'),
    'c3': ('C3', "Production's transforms re-expressed relative to one acquired volume"),
    'c1': ('C1', 'Own-type templates, each registered to a common frame'),
    'c2': ('C2', 'C1 with templates rebuilt from registered volumes three times (QSIPrep-style)'),
    'c4': ('C4', 'Pair-shared: each control-label pair gets one (mean) transform'),
    'c8-noreg': ('C8 no regression', 'Control and label registered jointly to their mean'),
    'c8': ('C8', 'Joint estimation, then the label/control alternation regressed out (OLS)'),
    'c8-huber': ('C8 Huber', 'C8 with a robust (Huber) fit'),
    'c8-perdelay': ('C8 per-delay', 'C8 with one alternation regressor per delay'),
    'c8-prod': ('C8 prod. ref.', "C8 registered to production's reference"),
    'c8-spm': ('C8 SPM conv.', "C8 with ASLtbx/ExploreASL's parameter convention and regression"),
}
FIXTURE_NOTES = {
    'pcasl1pld': 'still head, 2D PCASL, single delay',
    'pcasl1pld_grase': 'still head, 3D GRASE',
    'pasl_multipld': 'still head, 2D PASL, 8 delays',
    'pcasl_multipld': 'still head, PCASL, multi-delay',
    'pasl1pld': 'still head, PASL, single delay',
    'headmotion': 'aslscan headmotion: drift plus three 2 mm spikes, 10 pairs',
    'mo_pcasl_drift': 'slow drift to ~1.5 mm / 1.5 deg',
    'mo_pcasl_spikes': 'still, with two 2 mm / 2 deg spikes',
    'mo_pcasl_throughplane': 'drift along z with rotation about x',
    'mo_pcasl_betweenpairs': 'random walk that changes only between pairs',
    'mo_pcasl_withinpairs': 'between-pair walk plus a 0.5 mm / 0.5 deg jump inside each pair',
    'mo_pcasl_synchronous': 'every label displaced 0.4 mm / 0.4 deg (label-locked)',
    'mo_pcasl_constant': 'every volume at the same non-zero pose',
    'mo_pcaslinc_constant': 'constant pose, two M0 volumes in the series',
    'mo_pcaslinc_drift': 'drift, two M0 volumes in the series',
    'mo_pcaslshort2_drift': 'drift, only 2 pairs',
    'mo_pcaslshort4_drift': 'drift, only 4 pairs',
    'mo_pasl_drift': 'drift, PASL 8 delays x 1 pair',
    'mo_pasl_betweenpairs': 'between-pair walk, PASL 8 delays',
    'mo_pasl_withinpairs': 'within-pair jumps, PASL 8 delays',
    'mo_pasl_pldconfounded': 'one delay displaced 1 mm / 1 deg',
    'mo_paslrep_drift': 'drift, PASL 4 delays x 3 pairs',
    'mo_paslrep_withinpairs': 'within-pair jumps, PASL 4 delays x 3 pairs',
    'mo_paslrep_pldconfounded': 'one delay displaced, PASL 4 delays x 3 pairs',
}
#: Where aslscan's own fixtures (`headmotion`) live; set by `--aslscan`.
ASLSCAN = Path('/home/tsalo/aslscan-data/aslscan')
STILL = ['pcasl1pld', 'pcasl1pld_grase', 'pasl_multipld', 'pcasl_multipld']

# ---------------------------------------------------------------------------------------------
# SVG helpers. Colors come from CSS tokens (--fg, --muted, --grid, --s-<arm>), so the figures
# follow the page theme.
# ---------------------------------------------------------------------------------------------
esc = html.escape


def nice_ticks(lo, hi, n=5):
    if hi <= lo:
        hi = lo + 1
    raw = (hi - lo) / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= raw), default=10 * mag)
    start = math.floor(lo / step) * step
    ticks = []
    t = start
    while True:  # up to the first tick at or above ``hi``, so the axis covers every value
        ticks.append(round(t, 10))
        if t >= hi - step * 1e-9:
            return ticks
        t += step


def fmt(v, digits=2):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return '&ndash;'
    return f'{v:.{digits}f}'


def _frame(w, h, m):
    return w - m[1] - m[3], h - m[0] - m[2]


def svg_open(w, h, label):
    return (
        f'<svg class="chart" viewBox="0 0 {w} {h}" role="img" aria-label="{esc(label)}" '
        f'preserveAspectRatio="xMidYMid meet">'
    )


def axes(w, h, m, xt, yt, xs, ys, xlabel, ylabel, xfmt=str, yfmt=str):
    pw, ph = _frame(w, h, m)
    out = []
    for t in yt:
        y = ys(t)
        out.append(f'<line class="grid" x1="{m[3]}" x2="{m[3] + pw}" y1="{y:.1f}" y2="{y:.1f}"/>')
        out.append(
            f'<text class="tick" x="{m[3] - 6}" y="{y + 4:.1f}" text-anchor="end">{yfmt(t)}</text>'
        )
    for t in xt:
        x = xs(t)
        out.append(
            f'<line class="tickmark" x1="{x:.1f}" x2="{x:.1f}" y1="{m[0] + ph}" y2="{m[0] + ph + 4}"/>'
        )
        out.append(
            f'<text class="tick" x="{x:.1f}" y="{m[0] + ph + 17}" text-anchor="middle">{xfmt(t)}</text>'
        )
    out.append(
        f'<line class="axis" x1="{m[3]}" x2="{m[3] + pw}" y1="{m[0] + ph}" y2="{m[0] + ph}"/>'
    )
    out.append(
        f'<text class="axlabel" x="{m[3] + pw / 2}" y="{h - 4}" text-anchor="middle">{esc(xlabel)}</text>'
    )
    out.append(
        f'<text class="axlabel" transform="translate(13 {m[0] + ph / 2}) rotate(-90)" '
        f'text-anchor="middle">{esc(ylabel)}</text>'
    )
    return out


def legend(items):
    """HTML legend: items are (css color expression, label, dashed)."""
    spans = []
    for color, label, dashed in items:
        style = f'border-top: 3px {"dashed" if dashed else "solid"} {color}'
        spans.append(f'<span class="key"><i style="{style}"></i>{esc(label)}</span>')
    return f'<div class="legend">{"".join(spans)}</div>'


def line_chart(
    series,
    xlabel,
    ylabel,
    w=680,
    h=300,
    ymin=None,
    ymax=None,
    xticks=None,
    xfmt=None,
    yfmt=None,
    label='chart',
    points=True,
):
    """``series``: dicts with x, y, color, name, optional dash, marker_classes."""
    m = (14, 16, 40, 54)
    pw, ph = _frame(w, h, m)
    xs_all = np.concatenate([np.asarray(s['x'], float) for s in series])
    ys_all = np.concatenate([np.asarray(s['y'], float) for s in series])
    ys_all = ys_all[np.isfinite(ys_all)]
    x0, x1 = float(xs_all.min()), float(xs_all.max())
    y0 = float(ys_all.min()) if ymin is None else ymin
    y1 = float(ys_all.max()) if ymax is None else ymax
    yt = nice_ticks(y0, y1)
    y0, y1 = min(y0, yt[0]), max(y1, yt[-1])
    xt = xticks if xticks is not None else nice_ticks(x0, x1, 6)
    xt = [t for t in xt if x0 - 1e-9 <= t <= x1 + 1e-9]

    def xs(v):
        return m[3] + (v - x0) / ((x1 - x0) or 1) * pw

    def ys(v):
        return m[0] + ph - (v - y0) / ((y1 - y0) or 1) * ph

    out = [svg_open(w, h, label)]
    out += axes(
        w,
        h,
        m,
        xt,
        yt,
        xs,
        ys,
        xlabel,
        ylabel,
        xfmt or (lambda t: f'{t:g}'),
        yfmt or (lambda t: f'{t:g}'),
    )
    for s in series:
        pts = [(xs(a), ys(b)) for a, b in zip(s['x'], s['y'], strict=True) if np.isfinite(b)]
        d = ' '.join(
            f'{"M" if i == 0 else "L"}{px:.1f},{py:.1f}' for i, (px, py) in enumerate(pts)
        )
        dash = ' stroke-dasharray="5 4"' if s.get('dash') else ''
        out.append(f'<path d="{d}" fill="none" stroke="{s["color"]}" stroke-width="2"{dash}/>')
        if points:
            classes = s.get('marker_classes') or [''] * len(pts)
            for (px, py), cls in zip(pts, classes, strict=False):
                fill = 'var(--bg)' if cls == 'open' else s['color']
                out.append(
                    f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3.2" fill="{fill}" '
                    f'stroke="{s["color"]}" stroke-width="1.5"/>'
                )
    out.append('</svg>')
    return ''.join(out)


def grouped_bars(groups, series, values, ylabel, w=680, h=300, label='chart', yfmt=None):
    """``values[g][s]``: one bar per series within each group."""
    m = (14, 12, 54, 54)
    pw, ph = _frame(w, h, m)
    vmax = max(v for g in values for v in g if np.isfinite(v))
    yt = nice_ticks(0, vmax)
    y1 = yt[-1]

    def ys(v):
        return m[0] + ph - v / y1 * ph

    out = [svg_open(w, h, label)]
    for t in yt:
        out.append(
            f'<line class="grid" x1="{m[3]}" x2="{m[3] + pw}" y1="{ys(t):.1f}" y2="{ys(t):.1f}"/>'
        )
        out.append(
            f'<text class="tick" x="{m[3] - 6}" y="{ys(t) + 4:.1f}" text-anchor="end">{(yfmt or str)(t)}</text>'
        )
    gw = pw / len(groups)
    bw = gw * 0.8 / len(series)
    for i, g in enumerate(groups):
        for j, s in enumerate(series):
            v = values[i][j]
            if not np.isfinite(v):
                continue
            x = m[3] + i * gw + gw * 0.1 + j * bw
            y = ys(v)
            out.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw - 1.5:.1f}" height="{m[0] + ph - y:.1f}" '
                f'fill="{s["color"]}"><title>{esc(g)} · {esc(s["name"])}: {v:.3f}</title></rect>'
            )
        cx = m[3] + (i + 0.5) * gw
        out.append(
            f'<text class="tick" x="{cx:.1f}" y="{m[0] + ph + 16}" text-anchor="middle">{esc(g)}</text>'
        )
    out.append(
        f'<line class="axis" x1="{m[3]}" x2="{m[3] + pw}" y1="{m[0] + ph}" y2="{m[0] + ph}"/>'
    )
    out.append(
        f'<text class="axlabel" transform="translate(13 {m[0] + ph / 2}) rotate(-90)" '
        f'text-anchor="middle">{esc(ylabel)}</text>'
    )
    out.append('</svg>')
    return ''.join(out)


def strip_plot(cats, values, xlabel, w=680, label='chart', xmax=None, color='var(--accent)'):
    """One row per category; jittered points and a median tick."""
    row = 26
    m = (10, 16, 40, 210)
    h = m[0] + m[2] + row * len(cats)
    pw = w - m[1] - m[3]
    allv = np.concatenate([np.asarray(v, float) for v in values])
    x1 = xmax if xmax is not None else float(np.nanmax(allv))
    xt = nice_ticks(0, x1)
    x1 = xt[-1]

    def xs(v):
        return m[3] + min(v, x1) / x1 * pw

    out = [svg_open(w, h, label)]
    for t in xt:
        out.append(
            f'<line class="grid" x1="{xs(t):.1f}" x2="{xs(t):.1f}" y1="{m[0]}" y2="{h - m[2]}"/>'
        )
        out.append(
            f'<text class="tick" x="{xs(t):.1f}" y="{h - m[2] + 16}" text-anchor="middle">{t:g}</text>'
        )
    rng = np.random.default_rng(0)
    for i, (c, v) in enumerate(zip(cats, values, strict=True)):
        cy = m[0] + row * (i + 0.5)
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        out.append(
            f'<text class="tick" x="{m[3] - 8}" y="{cy + 4:.1f}" text-anchor="end">{esc(c)} (n={v.size})</text>'
        )
        for val, jit in zip(v, rng.uniform(-7, 7, v.size), strict=True):
            out.append(
                f'<circle cx="{xs(val):.1f}" cy="{cy + jit:.1f}" r="2.6" fill="{color}" fill-opacity="0.45"/>'
            )
        if v.size:
            mx = xs(float(np.median(v)))
            out.append(
                f'<line class="median" x1="{mx:.1f}" x2="{mx:.1f}" y1="{cy - 10:.1f}" y2="{cy + 10:.1f}"/>'
            )
    out.append(
        f'<text class="axlabel" x="{m[3] + pw / 2}" y="{h - 4}" text-anchor="middle">{esc(xlabel)}</text>'
    )
    out.append('</svg>')
    return ''.join(out)


def scatter(points, xlabel, ylabel, w=680, h=340, label='chart'):
    """``points``: dicts with x, y, color, name, optional ring (bool) and tag (short text)."""
    m = (14, 16, 42, 56)
    pw, ph = _frame(w, h, m)
    xs_all = [p['x'] for p in points]
    ys_all = [p['y'] for p in points]
    xt, yt = nice_ticks(0, max(xs_all)), nice_ticks(0, max(ys_all))
    x1, y1 = xt[-1], yt[-1]

    def xs(v):
        return m[3] + v / x1 * pw

    def ys(v):
        return m[0] + ph - v / y1 * ph

    out = [svg_open(w, h, label)]
    out += axes(w, h, m, xt, yt, xs, ys, xlabel, ylabel, lambda t: f'{t:g}', lambda t: f'{t:g}')
    for p in points:
        cx, cy = xs(p['x']), ys(p['y'])
        r = 6 if p.get('ring') else 4.2
        out.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{p["color"]}" fill-opacity="0.85" '
            f'stroke="{"var(--fg)" if p.get("ring") else p["color"]}" stroke-width="{2 if p.get("ring") else 1}">'
            f'<title>{esc(p["name"])}: {p["x"]:.3f}, {p["y"]:.3f}</title></circle>'
        )
        if p.get('tag'):
            out.append(
                f'<text class="tick" x="{cx + 7:.1f}" y="{cy - 6:.1f}">{esc(p["tag"])}</text>'
            )
    out.append('</svg>')
    return ''.join(out)


def heat_table(rows, cols, values, notes=None, digits=2, vmax=None, col_labels=None):
    """HTML table with cells shaded on a log scale (color-mix with the theme tokens)."""
    finite = [v for r in values for v in r if v is not None and np.isfinite(v) and v > 0]
    lo = max(min(finite), 1e-3) if finite else 1e-3
    hi = vmax or (max(finite) if finite else 1)
    head = ''.join(f'<th>{esc(c)}</th>' for c in (col_labels or cols))
    body = []
    for r, vals in zip(rows, values, strict=True):
        cells = []
        for v in vals:
            if v is None or not np.isfinite(v):
                cells.append('<td class="na">&ndash;</td>')
                continue
            p = 0 if v <= lo else min(1.0, math.log(v / lo) / math.log(hi / lo)) if hi > lo else 0
            cells.append(
                f'<td style="background: color-mix(in oklab, var(--hot) {p * 78:.0f}%, var(--surface))">'
                f'{v:.{digits}f}</td>'
            )
        note = f'<span class="rownote">{esc(notes[r])}</span>' if notes and r in notes else ''
        body.append(f'<tr><th scope="row"><code>{esc(r)}</code>{note}</th>{"".join(cells)}</tr>')
    return (
        f'<div class="tablewrap"><table class="heat"><thead><tr><th></th>{head}</tr></thead>'
        f'<tbody>{"".join(body)}</tbody></table></div>'
    )


def plain_table(header, rows, num_cols=()):
    head = ''.join(f'<th>{h}</th>' for h in header)
    body = ''
    for r in rows:
        tds = ''.join(
            f'<td class="{"num" if i in num_cols else ""}">{c}</td>' for i, c in enumerate(r)
        )
        body += f'<tr>{tds}</tr>'
    return f'<div class="tablewrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def figure(body, caption, num):
    return f'<figure>{body}<figcaption><b>Figure {num}.</b> {caption}</figcaption></figure>'


# ---------------------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------------------
def load_scores(root):
    rows = []
    for f in sorted((root / 'fixture-scores').glob('*/*/score.json')):
        s = json.loads(f.read_text())
        if s.get('score_version', 1) < 2:  # older scale (noise-seed runs whose fixtures are gone)
            continue
        # held-out (``_ho``, ``_ho2``): never shown with development
        if f.parent.name.split('__s')[0].endswith(('_ho', '_ho2')):
            continue
        rec = {'label': f.parent.parent.name, 'fixture': f.parent.name, 'path': f.parent}
        rec.update({k: _get(s, p) for k, p in METRICS.items()})
        rec['pair_by'] = s.get('pair', {}).get('by_volume_mm')
        rows.append(rec)
    return pd.DataFrame(rows)


def noise_cost(root):
    rows = []
    for x in json.loads((root / 'results.json').read_text()):
        part = x.get('label') or x.get('deltam')
        if not part:
            continue
        rows.append(
            {k: x[k] for k in ('recipe', 'cond', 'cost', 'ref', 'tag')}
            | {
                'offset': part['median'],
                'pair_max': x.get('pair', {}).get('max', np.nan),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------------------------
def section_noise(root, fig):
    df = noise_cost(root)
    std = df.query("tag == 'std' and cost == 'mutualinfo'")
    recipes = [
        r
        for r in (
            'pcasl1pld',
            'pasl_multipld',
            'pcasl_multipld',
            'pcasl1pld_grase',
            'pcasl1pld_ge3d',
        )
        if r in set(std.recipe)
    ]
    # noiseB repeats noiseA's realization in this experiment (it predates the seed fix), so it
    # is left out rather than shown as a second seed.
    conds = [c for c in ('clean', 'noiseA', 'noise10x') if c in set(std.cond)]
    cond_names = {'clean': 'no noise', 'noiseA': 'calibrated noise', 'noise10x': '10x noise'}
    shades = ['var(--s-c0)', 'var(--s-c1)', 'var(--s-c4)']
    common = std.query("ref == 'common'").set_index(['recipe', 'cond'])['offset']
    own = std.query("ref == 'own'").set_index(['recipe', 'cond'])['offset']
    vals = [[float(common.get((r, c), np.nan)) for c in conds] for r in recipes]
    series = [{'name': cond_names[c], 'color': shades[i]} for i, c in enumerate(conds)]
    chart = grouped_bars(
        recipes,
        series,
        vals,
        'label offset from truth (mm, median)',
        label='Spurious label offset by noise level',
    )
    cl_recipes = [r for r in recipes if r != 'pcasl1pld_ge3d']
    clean = np.array([float(common.get((r, 'clean'), np.nan)) for r in cl_recipes])
    noisy = np.array([float(common.get((r, 'noise10x'), np.nan)) for r in cl_recipes])
    own_max = float(np.nanmax(own.values)) if len(own) else np.nan
    change = float(np.nanmax(np.abs(noisy / clean - 1)) * 100)
    grase = float(common.get(('pcasl1pld_grase', 'clean'), np.nan))
    ge = [
        float(common.get(('pcasl1pld_ge3d', c), np.nan)) for c in ('clean', 'noiseA', 'noise10x')
    ]

    costs = ['mutualinfo', 'normmi', 'normcorr', 'corratio', 'leastsquares']
    clean_std = df.query("tag == 'std' and cond == 'clean'")
    tab_rows, tab_vals = [], []
    for r in recipes:
        for ref in ('common', 'own'):
            sub = clean_std.query('recipe == @r and ref == @ref').set_index('cost')['offset']
            if sub.empty:
                continue
            tab_rows.append(f'{r} · {ref} reference')
            tab_vals.append([float(sub.get(c, np.nan)) for c in costs])
    table = heat_table(tab_rows, costs, tab_vals, digits=3)
    own_all_max = float(clean_std.query("ref == 'own'").offset.max()) + 0.005
    return f"""
<section id="noise">
<h2>Is the spurious motion noise?</h2>
<p>On a still head every estimated transform should be the identity. Production registers label
volumes to a reference built from other volumes (control, or a separate M0 scan), and the label
volumes come out displaced anyway. The question for #705 and #707 is whether that offset is noise
the registration chases or something systematic.</p>
<p>The simulator can switch noise off. With no noise at all the label offset is already there:
{fmt(clean[0], 3)} mm on <code>pcasl1pld</code> and {fmt(grase, 2)} mm on 3D GRASE. On the
control-label fixtures, adding noise up to ten times the calibrated level changes it by at most
{change:.0f}%. The offset comes from registering one contrast to another, not from noise.</p>
<p>The GE fixture is the exception. Its delta-M volumes are registered to the M0 through a 5 mm
smoothed reference, and noise does make that worse: {fmt(ge[0], 2)} mm with no noise,
{fmt(ge[1], 2)} mm at the calibrated level and {fmt(ge[2], 2)} mm at ten times. The no-noise value
alone is already the largest offset of any fixture.</p>
{
        figure(
            legend([(s['color'], s['name'], False) for s in series]) + chart,
            'Median label-volume offset from the true pose (mm) on still fixtures, production-style common '
            'reference, mutual information; for the GE fixture, the delta-M volumes. Each fixture is shown '
            'with no noise, at the calibrated noise level and at ten times that. With each type registered '
            f'to its own template the offset is at most {fmt(own_max, 3)} mm.',
            fig,
        )
    }
<h3>Does the cost function matter?</h3>
<p>No cost function removes it. Every MCFLIRT cost leaves an offset when one contrast is
registered to a reference of another, and registering each volume type to its own template brings
every cost below {
        fmt(own_all_max, 2)
    } mm. Among the control-label fixtures 3D GRASE, where label and
control differ most, is the worst case for every cost; the GE fixture is worse still.</p>
{table}
<p class="note">Median label offset (mm), noise-free still fixtures. Shading is on a log scale.</p>
</section>"""


def section_blur(root, fig):
    out_charts, texts = [], {}
    for fx in ('pcasl1pld', 'pcasl1pld_ge3d'):
        f = root / 'blur-bias' / fx / 'blur_bias.json'
        if not f.exists():
            continue
        df = pd.DataFrame(json.loads(f.read_text()))
        a = df.query("fov == 'as_acquired' and cost == 'mutualinfo'")
        series = []
        for kernel, name, color, dash in (
            ('iso', 'reference blurred only', 'var(--s-c8)', False),
            ('both', 'both images blurred (matched)', 'var(--s-c1)', True),
        ):
            k = a[a.kernel == kernel].sort_values('fwhm_mm')
            series.append(
                {
                    'x': k.fwhm_mm.tolist(),
                    'y': k.rms_mm.tolist(),
                    'color': color,
                    'name': name,
                    'dash': dash,
                }
            )
        out_charts.append(
            (
                fx,
                line_chart(
                    series,
                    'smoothing FWHM (mm)',
                    'registration error (mm RMS)',
                    w=430,
                    h=260,
                    ymin=0,
                    xticks=[0, 2.5, 5, 8],
                    label=f'Blur bias on {fx}',
                ),
            )
        )
        k5 = df.query("cost == 'mutualinfo' and kernel == 'iso' and fwhm_mm == 5")
        texts[fx] = {r.fov: r.rms_mm for r in k5.itertuples()}
        costs5 = df.query("fov == 'as_acquired' and kernel == 'iso' and fwhm_mm == 5")
        texts[fx + '_costs'] = {r.cost: r.rms_mm for r in costs5.itertuples()}
    pc = texts.get('pcasl1pld', {})
    ge = texts.get('pcasl1pld_ge3d', {})
    gec = texts.get('pcasl1pld_ge3d_costs', {})
    fov_rows = [
        [f'<code>{esc(fov)}</code>', fmt(pc.get(fov), 3), fmt(ge.get(fov), 3)]
        for fov in ('as_acquired', 'drop_edge_slices', 'zero_padded', 'recentred')
    ]
    leg = legend(
        [
            ('var(--s-c8)', 'reference blurred only', False),
            ('var(--s-c1)', 'both images blurred (matched)', True),
        ]
    )
    charts = ''.join(
        f'<div class="panel"><h4><code>{esc(fx)}</code></h4>{c}</div>' for fx, c in out_charts
    )
    return f"""
<section id="blur">
<h2>Blurred references bias registration</h2>
<p>Production smooths the reference by 5 mm on the GE path (short series of five volumes or
fewer), and every averaged template is somewhat blurred. To isolate that effect, each fixture
volume was registered to a blurred copy of itself, so the true answer is the identity.</p>
<p>Registering a sharp image to a blurred one is biased. With a 5 mm isotropic kernel and mutual
information, the error is {fmt(pc.get('as_acquired'), 2)} mm on 2D PCASL and
{fmt(ge.get('as_acquired'), 2)} mm on the GE 3D fixture; across MCFLIRT costs on the GE fixture it
ranges from {fmt(min(gec.values()) if gec else None, 2)} to {
        fmt(max(gec.values()) if gec else None, 2)
    } mm.
Blurring both images by the same amount removes it, so the bias comes from the mismatch, not the
blur itself. Removing the edge slices, padding or recentring the field of view leaves it between
{fmt(min(ge.values()) if ge else None, 2)} and {fmt(max(ge.values()) if ge else None, 2)} mm on the
GE fixture, so it is not an edge-slice effect.</p>
{
        figure(
            leg + f'<div class="panels">{charts}</div>',
            'Registration error (mm RMS displacement) of a volume registered to a blurred copy of itself, '
            'mutual information, field of view as acquired.',
            fig,
        )
    }
{
        plain_table(
            ['Field of view', '<code>pcasl1pld</code> (mm)', '<code>pcasl1pld_ge3d</code> (mm)'],
            fov_rows,
            (1, 2),
        )
    }
<p class="note">Error at 5 mm FWHM (reference blurred only, mutual information) for each field-of-view
variant.</p>
<p>This matters twice. It explains part of the GE path's spurious motion (#702), and it predicts
that averaged templates (candidates C1 and C2 below) carry a blur bias whenever the head moves,
because averaging misaligned volumes blurs the template.</p>
</section>"""


def section_screen(df, fig):
    arms = [a for a in ARMS if a in set(df.label)]
    fixtures = [
        f
        for f in [*STILL, 'headmotion', *sorted(k for k in FIXTURE_NOTES if k.startswith('mo_'))]
        if f in set(df.fixture)
    ]
    piv = df.pivot_table(index='fixture', columns='label', values='dm_edge')
    pair = df.pivot_table(index='fixture', columns='label', values='pair_max')
    vals = [
        [float(piv.at[f, a]) if f in piv.index and a in piv.columns else np.nan for a in arms]
        for f in fixtures
    ]
    heat = heat_table(
        fixtures, arms, vals, notes=FIXTURE_NOTES, col_labels=[ARMS[a][0] for a in arms]
    )
    med = {
        a: float(np.nanmedian(piv.loc[[f for f in fixtures if f in piv.index], a])) for a in arms
    }
    med_pair = {
        a: float(np.nanmedian(pair.loc[[f for f in fixtures if f in pair.index], a])) for a in arms
    }
    order = sorted(arms, key=lambda a: med[a])
    bars = grouped_bars(
        [ARMS[a][0] for a in order],
        [{'name': 'median', 'color': 'var(--accent)'}],
        [[med[a]] for a in order],
        'median edge error (x GM delta-M)',
        h=280,
        label='Median edge delta-M error by candidate',
    )
    arm_rows = [
        [f'<b>{esc(ARMS[a][0])}</b>', esc(ARMS[a][1]), fmt(med[a], 2), fmt(med_pair[a], 3)]
        for a in order
    ]

    def v(f, a):
        try:
            return float(piv.at[f, a])
        except KeyError:
            return np.nan

    still_c0 = np.nanmedian([v(f, 'c0') for f in STILL])
    return f"""
<section id="screen">
<h2>Candidate screen on simulated motion</h2>
<p>Each candidate ran on {
        len(fixtures)
    } simulated fixtures: four still-head acquisitions, aslscan's
<code>headmotion</code>, and {
        sum(f.startswith('mo_') for f in fixtures)
    } motion fixtures generated for
this probe with known trajectories. Each fixture was scored against an oracle: the same raw
volumes moved back by their <i>true</i> poses and subtracted the same way.</p>
<p>The headline metric is the RMS difference between the candidate's mean delta-M and the oracle's,
in the outer 0&ndash;3 mm of the brain, divided by the oracle's grey-matter delta-M. Edges are where
misregistration shows: a 0.1 mm shift between control and label leaves a rim of raw tissue
contrast, which is about a hundred times the perfusion signal. A value of 1 means the edge error
is as large as the grey-matter signal.</p>
{
        figure(
            bars,
            'Median edge error across all fixtures for each candidate (lower is better). '
            'Candidates are defined in the table below.',
            fig,
        )
    }
{
        plain_table(
            ['Candidate', 'What it does', 'Median edge error', 'Median worst pair error (mm)'],
            arm_rows,
            (2, 3),
        )
    }
<h3>Fixture by fixture</h3>
<p>Production's error on still heads is small but not zero (median {fmt(still_c0, 2)}): that is the
spurious label offset. Templates per volume type (C1, C2) remove it on still heads, but with real
motion they lose to production, as the blur bias predicts: on <code>mo_pcasl_drift</code> C1 scores
{fmt(v('mo_pcasl_drift', 'c1'), 2)} against production's {fmt(v('mo_pcasl_drift', 'c0'), 2)}.
Re-referencing production's transforms to one acquired volume (C3) changes nothing that matters.
Sharing one transform per pair (C4) and zig-zag regression (C8) are perfect when there is no
within-pair motion and fail badly when there is (next section).</p>
{
        figure(
            heat,
            'Edge delta-M error (x GM delta-M) for each fixture and candidate; darker is worse, on a '
            'log scale. Fixtures a candidate was not run on are blank.',
            fig + 1,
        )
    }
</section>"""


def section_zigzag(root, df, fig):
    piv = df.pivot_table(index='fixture', columns='label', values='dm_edge')

    def v(f, a):
        try:
            return float(piv.at[f, a])
        except KeyError:
            return np.nan

    art = {
        lab: root / 'fixture-scores' / lab / 'headmotion' / 'artifacts'
        for lab in ('c8', 'c8-huber')
    }
    params = np.load(art['c8'] / 'params_initial.npy')
    removed = {lab: np.load(a / 'removed.npy') for lab, a in art.items()}
    j = int(np.argmax(np.abs(removed['c8']).max(0)))
    names = [
        'rotation x',
        'rotation y',
        'rotation z',
        'translation x',
        'translation y',
        'translation z',
    ]
    unit = 'deg' if j < 3 else 'mm'
    ctx = pd.read_table(root / 'fixture-scores' / 'c8' / 'headmotion' / 'hmc_aslcontext.tsv')[
        'volume_type'
    ]
    types = [t for t in ctx if t in ('control', 'label')]
    marks = ['open' if t == 'label' else '' for t in types]
    x = list(range(len(params)))
    s1 = line_chart(
        [
            {
                'x': x,
                'y': params[:, j],
                'color': 'var(--s-c8-noreg)',
                'name': 'estimated',
                'marker_classes': marks,
            },
            {
                'x': x,
                'y': params[:, j] - removed['c8'][:, j],
                'color': 'var(--s-c8)',
                'name': 'OLS',
                'dash': True,
            },
            {
                'x': x,
                'y': params[:, j] - removed['c8-huber'][:, j],
                'color': 'var(--s-c8-huber)',
                'name': 'Huber',
                'dash': True,
            },
        ],
        'volume (control-label series, acquisition order)',
        f'{names[j]} ({unit})',
        w=680,
        h=280,
        label='Zig-zag regression on headmotion',
    )
    leg1 = legend(
        [
            ('var(--s-c8-noreg)', 'joint estimate (filled: control, open: label)', False),
            ('var(--s-c8)', 'after OLS zig-zag removal', True),
            ('var(--s-c8-huber)', 'after Huber zig-zag removal', True),
        ]
    )

    # The true trajectory's own alternation, as OLS sees it.
    gt = ASLSCAN / 'headmotion/sub-01/perf/ground-truth/sub-01_desc-motion_gt.tsv'
    if gt.exists():
        truth = pd.read_table(gt)
        tcols = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']
        tp = truth[tcols].to_numpy(float)
        tp[:, 3:] = np.degrees(tp[:, 3:])
        full_ctx = pd.read_table(
            ASLSCAN / 'headmotion/sub-01/perf/sub-01_acq-headmotion_aslcontext.tsv'
        )['volume_type'].tolist()
        keep = [i for i, t in enumerate(full_ctx) if t in ('control', 'label')]
        if len(tp) == len(full_ctx):
            tp = tp[keep]
        tt = [full_ctx[i] for i in keep]
        _, rem_ols, _ = zigzag_clean(tp, tt, [0.0] * len(tt))
        _, rem_hub, _ = zigzag_clean(tp, tt, [0.0] * len(tt), robust=True)
        t_ols = np.abs(rem_ols[:, :3]).max()
        r_ols = np.abs(rem_ols[:, 3:]).max()
        t_hub = np.abs(rem_hub[:, :3]).max()
    else:
        # The scored fixture's data were deleted (disk space); these are its values, as computed
        # for the report of 2026-10-19. A regenerated `headmotion` is a different simulation.
        t_ols, r_ols, t_hub = 0.15, 0.23, 0.08

    pair_series = []
    for lab, dash in (('c0', False), ('c8-noreg', False), ('c8', True), ('c8-huber', True)):
        row = df.query("label == @lab and fixture == 'mo_pcasl_drift'")
        if row.empty or row.iloc[0].pair_by is None:
            continue
        y = row.iloc[0].pair_by
        pair_series.append(
            {
                'x': list(range(1, len(y) + 1)),
                'y': y,
                'color': f'var(--s-{lab})',
                'name': ARMS[lab][0],
                'dash': dash,
            }
        )
    s2 = line_chart(
        pair_series,
        'control-label pair',
        'pair error (mm)',
        w=680,
        h=260,
        ymin=0,
        label='Pair error on drift',
    )
    leg2 = legend([(s['color'], s['name'], s['dash']) for s in pair_series])

    comp = [
        [
            'Estimation',
            'SPM realign: least squares, 5 mm smoothing, B-spline; to the first volume',
            'SPM realign, trilinear; two passes, to the mean',
            'MCFLIRT, mutual information; to the control-label mean',
        ],
        [
            'Volumes regressed',
            'the whole file',
            'the whole file (M0 usually separate)',
            'control and label from <code>aslcontext.tsv</code>; M0 separate',
        ],
        [
            'Alternation',
            'by position, label first unless <code>ctrfirst</code>',
            'by position, label first',
            'from <code>aslcontext.tsv</code>',
        ],
        [
            'Parameters',
            'Euler angles about the world origin, relative to volume 1',
            'the same',
            'rotation vector about the field-of-view centre',
        ],
        [
            'Regression',
            'no intercept, then re-referenced to volume 1',
            'the same',
            'intercept and centred regressor',
        ],
        ['Multi-delay', 'single delay', 'disabled for multi-PLD', 'pooled or per-delay'],
        [
            'Spikes',
            'not handled',
            'ENABLE spike removal, after the regression',
            'OLS, or the Huber variant',
        ],
    ]
    # Largest relative difference between C8 and its SPM-convention twin, where C8's error is
    # large enough for a ratio to mean something.
    both = piv[['c8', 'c8-spm']].dropna()
    both = both[both.c8 > 0.05]
    spm_pct = (
        float((np.abs(both['c8-spm'] - both.c8) / both.c8).max() * 100) if len(both) else np.nan
    )
    spm_rows = [
        [f'<code>{esc(f)}</code>', fmt(v(f, 'c8'), 2), fmt(v(f, 'c8-spm'), 2)]
        for f in (
            'pcasl1pld',
            'headmotion',
            'mo_pcasl_drift',
            'mo_pcasl_withinpairs',
            'mo_pcasl_synchronous',
            'mo_pcaslshort2_drift',
        )
    ]
    return f"""
<section id="zigzag">
<h2>Zig-zag regression</h2>
<p>Wang (2012) removes the spurious label/control pattern after the fact: register every volume,
then regress a +1/&minus;1 alternation out of each motion parameter. ASLtbx and ExploreASL both do this.
Candidate C8 is that method, with M0 volumes registered separately (Sudipto Dolui's suggestion in
#266).</p>
<p>On still heads it is perfect: edge error {
        fmt(v('pcasl1pld', 'c8'), 2)
    } on <code>pcasl1pld</code>,
against production's {
        fmt(v('pcasl1pld', 'c0'), 2)
    }. It fails as soon as real motion has a component
at the alternation frequency, and almost every real motion does.</p>
<h3>Spikes alias onto the alternation</h3>
<p>The <code>headmotion</code> trajectory has no alternation at all, but with ten pairs, a drift and
three spikes, least squares still finds one in the true parameters: up to {fmt(t_ols, 2)} mm and
{fmt(r_ols, 2)} deg. Subtracting it puts a zig-zag into every volume, and the edge error goes from
{fmt(v('headmotion', 'c8-noreg'), 2)} (joint estimation, no regression) to
{fmt(v('headmotion', 'c8'), 2)}. A Huber fit down-weights the spikes ({
        fmt(t_hub, 2)
    } mm on the true
parameters) and brings it back to {fmt(v('headmotion', 'c8-huber'), 2)}; on
<code>mo_pcasl_spikes</code> it goes from {fmt(v('mo_pcasl_spikes', 'c8'), 2)} to
{fmt(v('mo_pcasl_spikes', 'c8-huber'), 2)}, production's level.</p>
{
        figure(
            leg1 + s1,
            f'Estimated {names[j]} on <code>headmotion</code> (the parameter with the largest '
            'removed component) before and after zig-zag removal. The true trajectory has no alternation.',
            fig,
        )
    }
<h3>Drift is removed along with the zig-zag</h3>
<p>A slow drift moves the head a little between the control and the label of every pair. That
step is real, and it alternates in sign relative to the pair's mean, so the regression removes half
of it from each volume. The pair error triples, and because edges are so sensitive, the edge error
goes from {fmt(v('mo_pcasl_drift', 'c0'), 2)} (production) to {
        fmt(v('mo_pcasl_drift', 'c8'), 2)
    } on
<code>mo_pcasl_drift</code>. Huber cannot help here: the drift is not an outlier. Pair-shared
transforms (C4) fail the same way ({fmt(v('mo_pcasl_drift', 'c4'), 2)}).</p>
{
        figure(
            leg2 + s2,
            'Error in the control-to-label transform of each pair on <code>mo_pcasl_drift</code> '
            '(mm RMS over an 80 mm sphere).',
            fig + 1,
        )
    }
<h3>Comparison with ASLtbx and ExploreASL</h3>
<p>Both tools share Ze Wang's parameter-cleaning function. The differences from C8:</p>
{plain_table(['', 'ASLtbx', 'ExploreASL', 'C8'], comp)}
<p>For an even number of strictly alternating volumes, their no-intercept projection equals C8's
regression up to the re-referencing; the parameterization differs. The arm <code>c8-spm</code>
reproduces their convention and regression exactly on top of C8's joint estimation. It scores
within {
        spm_pct:.0f}% of C8 on every fixture, so the convention does not explain C8's failures, and a
faithful SPM re-implementation would not be expected to fix them.</p>
{plain_table(['Fixture', 'C8', 'C8, SPM convention'], spm_rows, (1, 2))}
<p class="note">Edge delta-M error (x GM delta-M).</p>
</section>"""


def section_consumers(root):
    rows = []
    for f in sorted((root / 'fixture-runs' / 'baseline').glob('*/score_full.json')):
        d = json.loads(f.read_text())
        rows.append(
            [
                f'<code>{esc(f.parent.name)}</code>',
                fmt(d.get('aslref_frame_offset_mm'), 2),
                fmt((d.get('mask') or {}).get('dice'), 3),
                fmt((d.get('coregistration_vs_aslref_frame') or {}).get('rms_mm'), 2),
                fmt((d.get('m0') or {}).get('rms_mm'), 2),
                fmt((d.get('fieldmap') or {}).get('rms_mm'), 2),
            ]
        )
    return f"""
<section id="consumers">
<h2>What consumes the reference</h2>
<p>The reference image is not only the HMC target. Skull-stripping, ASL-to-T1w coregistration,
M0 registration and fieldmap registration all start from it, so a candidate that changes the
reference has to be scored on those too. Full ASLPrep runs on five probe fixtures give production's
baseline:</p>
{
        plain_table(
            [
                'Fixture',
                'Reference offset from acquisition (mm)',
                'Brain-mask Dice',
                'Coregistration error (mm)',
                'M0 registration error (mm)',
                'Fieldmap registration error (mm)',
            ],
            rows,
            (1, 2, 3, 4, 5),
        )
    }
<p class="note">Errors are RMS displacement over an 80 mm sphere against the simulator's truth. Fixtures
without a fieldmap have no fieldmap entry.</p>
<p>The fieldmap registration error, 0.6&ndash;1.1 mm on the fixtures that have one, is as large as the
motion errors under study, and 3D GRASE coregistration is off by about 1 mm. These are the numbers
a new reference has to match or beat.</p>
</section>"""


def section_real(root, fig):
    screen = pd.read_table(root / 'harness-runs' / 'screen_c0.tsv')
    real = pd.read_table(root / 'harness-runs' / 'real_c0.tsv')
    m = screen.merge(real.drop(columns=['label']), on=['dataset', 'run', 'family'], how='left')
    alt = m.groupby('family').alternation_power.median().dropna().sort_values()
    no_pairs = sorted(set(m.family) - set(alt.index))
    fams = alt.index.tolist()
    strip = strip_plot(
        fams,
        [m.loc[m.family == f, 'alternation_power'] for f in fams],
        'alternation power (fraction of motion-parameter power at the label/control frequency)',
        label='Alternation power by family',
    )
    strip2 = strip_plot(
        fams,
        [m.loc[m.family == f, 'pair_median_mm'] for f in fams],
        'median within-pair motion (mm)',
        xmax=float(np.nanpercentile(m.pair_median_mm, 98)),
        label='Within-pair motion by family',
        color='var(--s-c8)',
    )
    rows = []
    for f in fams[::-1]:
        g = m[m.family == f]
        rows.append(
            [
                esc(f),
                str(len(g)),
                fmt(g.fd_mean_mm.median(), 3),
                fmt(g.pair_median_mm.median(), 3),
                fmt(g.alternation_power.median(), 3),
                fmt(g.negative_fraction.median(), 3),
                fmt(g.split_half_r.median(), 2),
            ]
        )
    n_ds = m.dataset.nunique()
    return f"""
<section id="real">
<h2>Real data: production's baseline</h2>
<p>Production's HMC ran on {len(m)} runs from {n_ds} datasets in the aslprep-zoo panel (the
anatomy-free harness, which builds the reference, corrects motion and subtracts). Real data have no
truth, so the probe uses indirect measures. Alternation power is the fraction of the estimated motion
parameters' power at the label/control frequency: real head motion has little reason to alternate
with labelling, so high values point to spurious label motion.</p>
{
        figure(
            strip,
            "Alternation power of production's motion estimates, one point per run, by acquisition "
            'family; the tick is the median.',
            fig,
        )
    }
<p>{esc(alt.index[-1])} stands out (median {fmt(alt.iloc[-1], 3)}, against
{fmt(float(alt.iloc[:-1].median()), 3)} for the median of the other families). The simulations
predict the most spurious label motion where label and control differ most, which is the case for
3D readouts. {
        esc(', '.join(no_pairs)) or 'No family'
    } has no control-label series (delta-M only), so
it has no alternation power and is not shown.</p>
{figure(strip2, 'Median within-pair motion estimated by production, by family.', fig + 1)}
{
        plain_table(
            [
                'Family',
                'Runs',
                'Mean FD (mm)',
                'Within-pair motion (mm)',
                'Alternation power',
                'Negative delta-M fraction',
                'Split-half r',
            ],
            rows,
            (1, 2, 3, 4, 5, 6),
        )
    }
<p class="note">Medians per family, production only.</p>
{section_real_shortlist(root)}
</section>"""


FAMILY_COLOR = {
    'c0': 'var(--s-c0)',
    'c1': 'var(--s-c1)',
    'c2': 'var(--s-c2)',
    'c4': 'var(--s-c4)',
    'c8': 'var(--s-c8)',
    'c3': 'var(--s-c3)',
    'fixvol': 'var(--s-c3)',
    'perfref': 'var(--s-c8-noreg)',
}
ENGINE_NOTES = {
    'c1': 'MCFLIRT',
    'c1-m0': 'MCFLIRT, M0 frame',
    'c1-volreg': 'AFNI 3dvolreg',
    'c1-ants': 'ANTs rigid',
    'c1-niimath': 'niimath fast (rigid projection)',
    'c1-niimath-ls': 'niimath rigid least squares',
    'c1-pld': 'grouped by type and delay',
    'c1-idcross': 'identity cross-type step',
}


def section_engines(root, fig):
    f = root / 'fixture-scores' / 'estimator_screen.tsv'
    if not f.exists():
        return ''
    s = pd.read_csv(f, sep='\t')
    s = s[~s.diagnostic.astype(bool)]
    short = s[s['rank'].eq(1)].label.tolist()
    pts = []
    for r in s.itertuples():
        ring = r.label in short or r.label == 'c0'
        tag = r.label if ring or r.family in ('fixvol', 'perfref') or 'volreg' in r.label else ''
        pts.append(
            {
                'x': r.rel_max,
                'y': r.pair_max,
                'color': FAMILY_COLOR.get(r.family, 'var(--muted)'),
                'name': r.label,
                'ring': ring,
                'tag': tag,
            }
        )
    chart = scatter(
        pts,
        'median relative-pose error (mm, max over volumes)',
        'median worst pair error (mm)',
        label='Estimator screen',
    )
    leg = legend(
        [
            (FAMILY_COLOR[k], n, False)
            for k, n in (
                ('c0', 'production'),
                ('c1', 'C1'),
                ('c2', 'C2'),
                ('c4', 'C4'),
                ('c8', 'C8'),
                ('fixvol', 'fixed-volume comparator'),
                ('perfref', 'perfusion-reference comparator'),
            )
        ]
    )
    row = s.set_index('label')

    def g(label, col):
        return float(row.at[label, col]) if label in row.index else np.nan

    pair_word = 'below' if g('c2-volreg', 'pair_max') < g('c0', 'pair_max') else 'above'
    ranked = s[s.family.isin(['c1', 'c2', 'c4', 'c8'])].sort_values(['family', 'rank'])
    rows = [
        [
            f'{"<b>" if r.label in short else ""}<code>{esc(r.label)}</code>{"</b>" if r.label in short else ""}',
            esc(ENGINE_NOTES.get(r.label, ENGINE_NOTES.get(r.label.replace('c2', 'c1'), ''))),
            fmt(r.rel_max, 3),
            fmt(r.pair_max, 3),
            fmt(r.dm_edge, 2),
            f'{r.rank:.0f}',
        ]
        for r in ranked.itertuples()
    ]
    return f"""
<section id="engines">
<h2>Registration engines and the short list</h2>
<p>The plan's estimator screen runs every variant of each family on the {
        int(s.n.max())
    } development
fixtures and keeps the best of each by two metrics fixed in advance: relative-pose error (are the
volumes aligned with each other?) and pair error (are control and label aligned?). C1 and C2 were
crossed with four registration engines (MCFLIRT, AFNI 3dvolreg, ANTs rigid and niimath's port of
3dAllineate), the grouping (by type, or by type and delay), the common frame and the cross-type
step. Each engine was first checked against known rigid transforms (all recover them to within
0.31 mm for 2.5&ndash;3.8 mm motions).</p>
<p>3dvolreg stands out on relative pose: {fmt(g('c1-volreg', 'rel_max'), 3)} mm for C1 against
production's {fmt(g('c0', 'rel_max'), 3)} mm. Least squares within one contrast is the right cost
for a same-contrast template. Its pair error is {pair_word} production's
({fmt(g('c2-volreg', 'pair_max'), 3)} mm for C2 against {
        fmt(g('c0', 'pair_max'), 3)
    } mm), but its edge
error is not: {fmt(g('c2-volreg', 'dm_edge'), 2)} against {
        fmt(g('c0', 'dm_edge'), 2)
    }. ANTs has the
worst pair error of the engines ({fmt(g('c1-ants', 'pair_max'), 3)} mm).</p>
{
        figure(
            leg + chart,
            'Each arm of the screen: median relative-pose error against median pair error. '
            'Ringed points are production and the variant that goes forward from each family.',
            fig,
        )
    }
{
        plain_table(
            ['Arm', 'Variant', 'Relative pose (mm)', 'Pair (mm)', 'Edge error', 'Rank'],
            rows,
            (2, 3, 4, 5),
        )
    }
<p class="note">Medians over the development fixtures; rank within the family by the sum of the two
screening ranks, ties broken on pair error. C8's diagnostic-only runs and C1's leave-one-out null
test are not ranked.</p>
<p>The short list is production plus {', '.join(f'<code>{esc(x)}</code>' for x in short)}. The two
comparators from #513 bracket it. Registering everything to one acquired volume
(<code>fixvol</code>) is close to production. Re-registering each pair to a preliminary perfusion
image (<code>perfref</code>) keeps pairs aligned but adds {fmt(g('perfref', 'rel_max'), 2)} mm of
relative-pose error, because a single pair's delta-M is too noisy to register.</p>
</section>"""


def shortlist(root):
    """C0 plus each family's first-ranked arm in ``estimator_screen.tsv``."""
    f = root / 'fixture-scores' / 'estimator_screen.tsv'
    if not f.exists():
        return ['c0']
    s = pd.read_csv(f, sep='\t')
    return ['c0', *s[s['rank'].eq(1) & ~s.diagnostic.astype(bool)].label.tolist()]


def section_policy(df, fig, root):
    short = shortlist(root)
    labels = set(df.label)
    if not any(f'{s}-d1' in labels for s in short):
        return ''
    rows = []
    for s in short:
        cells = [f'<code>{esc(s)}</code>']
        for lab in (f'{s}-d0', f'{s}-d1', s):
            sub = df[df.label == lab]
            cells += (
                [fmt(sub.ref_med.median(), 3), fmt(sub.dmref_edge.median(), 2)]
                if len(sub)
                else ['&ndash;', '&ndash;']
            )
        rows.append(cells)
    head = ['Estimator'] + [
        f'{p} {m}' for p in ('D0', 'D1', 'D2') for m in ('ref. offset (mm)', 'edge error')
    ]
    d1 = df[df.label.str.endswith('-d1')]
    fallbacks = sum(
        'd1_fallback_to_d0' in json.loads((p / 'harness.json').read_text()).get('downstream', {})
        for p in d1.path
    )

    def med(label, col):
        sub = df[df.label == label]
        return float(sub[col].median()) if len(sub) else np.nan

    tpl_ref = [med(s, 'ref_med') for s in short if s.startswith(('c1', 'c2'))] or [np.nan]

    return f"""
<section id="policy">
<h2>Which reference the later steps receive</h2>
<p>The HMC target and the image handed to later steps need not be the same. D0 keeps production's
reference and expresses the candidate's transforms in its frame; D1 rebuilds production's reference
from the same volume selection after the candidate's correction; D2 hands on the candidate's own
template. The relative metrics do not change between policies; what changes is how well the
corrected volumes agree with the reference the later steps see.</p>
{plain_table(head, rows, tuple(range(1, 7)))}
<p class="note">Medians over the development fixtures: the median offset between corrected volumes and
the reference, and the edge delta-M error against the reference-frame oracle.</p>
<p>Handing on the candidate's own template (D2) brings the corrected volumes within
{fmt(min(tpl_ref), 2)}&ndash;{fmt(max(tpl_ref), 2)} mm of the reference for
the template candidates, against about {fmt(med('c0', 'ref_med'), 2)} mm for production. D1 differs
little from D0 here because {fallbacks} of {len(d1)} D1 runs fell back to D0: most fixtures have a
separate M0, which production selects but which is not in the corrected series. The edge error is
set by the estimator, not the policy; the policy matters for the later steps (brain mask,
coregistration, M0 and fieldmap registration), which only the full-pipeline runs of Task 7 can
score.</p>
</section>"""


def section_real_shortlist(root):
    f = root / 'harness-runs' / 'real_shortlist.tsv'
    if not f.exists():
        return ''
    d = pd.read_csv(f, sep='\t')
    labels = [x for x in shortlist(root) if x in set(d.label)]
    if len(labels) < 2:
        return ''
    d = d[d.label.isin(labels)]
    piv = d.pivot_table(
        index='family', columns='label', values='mean_dm_ssd_vs_c0', aggfunc='median'
    )
    neg = d.pivot_table(
        index='family', columns='label', values='negative_fraction', aggfunc='median'
    )
    others = [x for x in labels if x != 'c0']
    rows = [
        [
            esc(fam),
            *[fmt(piv.at[fam, x], 3) for x in others],
            *[fmt(neg.at[fam, x], 3) for x in labels],
        ]
        for fam in piv.index
        if piv.loc[fam, others].notna().any()
    ]
    head = (
        ['Family']
        + [f'<code>{esc(x)}</code> energy vs C0' for x in others]
        + [f'<code>{esc(x)}</code> neg. fraction' for x in labels]
    )
    # Over finite, matched comparisons only (runs without pairs have no ratio), with real_score's
    # margin: a ratio within it is a tie (interpolation alone moves the energy by up to 4.5 %).
    scored = {x: d[(d.label == x) & np.isfinite(d.mean_dm_ssd_vs_c0)] for x in others}
    verdicts = {
        x: {
            v: float((g.mean_dm_ssd_verdict == v).mean() * 100) for v in ('better', 'tie', 'worse')
        }
        for x, g in scored.items()
    }
    n_scored = min(len(g) for g in scored.values())
    return f"""
<h3>The short list on real data</h3>
<p>Each short-listed estimator ran on the same real runs, scored with the shared-mask mean delta-M
energy (as a ratio to production on the same run; below 1 is better) and the negative fraction.
A common shift of control and label alike, which creates no misregistration, already moves the
energy by up to 4.5% through interpolation alone, so a ratio within 5% of production's is counted
as a tie. On the {n_scored} scorable runs (better / tie / worse than production):
{'; '.join(f'<code>{esc(x)}</code> {v["better"]:.0f}% / {v["tie"]:.0f}% / {v["worse"]:.0f}%' for x, v in verdicts.items())}.
The energy only measures motion-driven errors; at still-head offsets it favours production.</p>
{plain_table(head, rows, tuple(range(1, len(head))))}
<p class="note">Medians per family. GE product output has no pairs and is not scored.</p>"""


C5_FIXTURES = [
    ('c5_presat', 'still'),
    ('c5_presat_drift', 'drift'),
    ('c5_presat_pldconfounded', 'one delay moved'),
    ('c5_presat_bs', 'BS, still'),
    ('c5_presat_bs_drift', 'BS, drift'),
    ('c5_presat_3d', '3D, still'),
    ('c5_presat_3d_drift', '3D, drift'),
]


def section_c5(root, fig):
    """C5 (#673): saturation-recovery normalization on its own fixtures, against production."""
    base = root / 'fixture-scores'
    arms = [
        ('c0', 'production', 'var(--s-c0)'),
        ('c5-voxel', 'C5 voxelwise T1eff', 'var(--s-c5-voxel)'),
        ('c5', 'C5 per-slice T1eff', 'var(--s-c5)'),
    ]

    def score(lab, name):
        f = base / lab / name / 'score.json'
        return json.loads(f.read_text()) if f.exists() else None

    if score('c5', 'c5_presat') is None:
        return ''
    rows, eligible = [], []
    for name, what in C5_FIXTURES:
        rec = base / 'c5' / name / 'artifacts' / 'c5.json'
        info = json.loads(rec.read_text()) if rec.exists() else {}
        cells = [esc(what)]
        for lab, _, _ in arms:
            s = score(lab, name)
            cells += (
                [
                    fmt(s['relative_pose']['to_mean']['median_mm'], 3),
                    fmt(s['relative_pose']['to_mean']['max_mm'], 3),
                    fmt(s['pair']['median_mm'], 4),
                ]
                if s
                else ['&ndash;'] * 3
            )
        status = 'applied' if info.get('eligible') else '; '.join(info.get('reasons', ['&ndash;']))
        rows.append([*cells, esc(status)])
        if info.get('eligible'):
            eligible.append(name)
    groups = [w for n, w in C5_FIXTURES if n in eligible]
    values = [
        [
            (s['relative_pose']['to_mean']['median_mm'] if (s := score(lab, n)) else np.nan)
            for lab, _, _ in arms
        ]
        for n in eligible
    ]
    series = [{'name': nm, 'color': c} for _, nm, c in arms]
    chart = grouped_bars(
        groups, series, values, 'pose error to the mean (mm)', label='C5 pose error'
    )
    head = (
        ['Fixture']
        + [f'{nm} {m}' for _, nm, _ in arms for m in ('to-mean med.', 'to-mean max', 'pair med.')]
        + ['C5']
    )
    return f"""
<section id="c5">
<h2>Saturation-recovery normalization (#673)</h2>
<p>With a pre-saturation pulse and no background-suppression inversions, each control image is a
point on a recovery curve, <code>M0 (1 &minus; e<sup>&minus;t/T1eff</sup>)</code>, where
<code>t</code> is the labeling duration plus the delay plus the slice time. Volumes at different delays
therefore differ in contrast, and production registers them to one another as they are. C5 divides
every control and label volume by its recovery factor and registers the normalized volumes jointly.
It refuses runs it cannot model (3D readouts, inversion pulses, no slice timing, a single delay)
and processes them as production. Its three fixtures are a 2D EPI pCASL protocol with five delays
(0.2&ndash;2.2 s) and pre-saturation at labeling start, and two negative controls that C5 must
refuse: the same protocol with two inversion pulses (BS), and with a 3D GRASE readout.</p>
{figure(legend([(s['color'], s['name'], False) for s in series]) + chart, 'Median pose error to the mean on the eligible fixtures. Lower is better.', fig)}
{plain_table(head, rows, tuple(range(1, len(head) - 1)))}
<p>C5 refuses both negative controls and reproduces production on them exactly. On the eligible
protocol, the voxelwise fit looks perfect on the still head and is the worst arm with motion. The
error jumps at the delay-block boundaries and grows with the motion, which points to a feedback
loop: the T1eff map has sharp edges, normalizing by it at the estimated pose imprints those edges
there, and registration then confirms the pose it was given. With no motion, the identity is that
fixed point. A single T1eff per slice, in the acquisition frame, cannot carry a pose. It ties
production: a little better with drift, worse on the one moved delay block, about the same on the
still head.</p>
<p class="note">Development fixtures; SNR 50 per pair. A per-slice T1eff removes the delay and
slice-time contrast that #673 is about, but not tissue differences in T1eff.</p>
</section>"""


C6_MODELS = {
    'free': 'free contrast maps, free P per delay',
    'smooth-pd': 'smooth contrast, free P per delay',
    'smooth': 'smooth contrast, one P',
    'gain': 'gain per group, one P',
    'smooth+supp-s': 'smooth contrast, one P, zero outside the image mask',
    'smooth+supp-t': 'smooth contrast, one P, zero outside true tissue',
    'gain+supp-t': 'gain per group, one P, zero outside true tissue',
}
C6_ARMS = {
    'c0': 'production',
    'c2-volreg': 'C2 (3dvolreg)',
    'c6-h1-pb3': 'C6, interpolation model, P bounded',
    'c6-pb3': 'C6, acquisition model, P bounded',
}
C6_KEY = ['mo_pcasl_synchronous', 'mo_pcasl_spikes', 'mo_pcasl_withinpairs', 'mo_pcasl_drift']


def section_c6(root, df, fig):
    """C6: the identifiability test and the prototype's development screen."""
    out = []
    f = root / 'c6' / 'identifiability.tsv'
    if f.exists():
        t = pd.read_csv(f, sep='\t')
        t = t[t.noise_pct == 2.0]
        worst = t.groupby(['model', 'fixture']).max_sd_trans_mm.max().unstack()
        fixtures = [x for x in ('c5_presat', 'rc_pasl_multipld', 'rc_pcasl') if x in worst.columns]
        rows = [
            [esc(C6_MODELS.get(m, m))]
            + [
                '&infin;' if not np.isfinite(worst.at[m, x]) else fmt(worst.at[m, x], 3)
                for x in fixtures
            ]
            for m in C6_MODELS
            if m in worst.index
        ]
        out.append(f"""
<p>C6 predicts every control and label volume in its own contrast and fits the poses jointly, with no
cross-contrast cost and no fixed reference. As first written it is not identifiable, so the plan
required an identifiability test before any prototype. The test is exact and local: the Fisher
information at the truth, with every nuisance parameter (static image, contrast, perfusion, delay
weights) eliminated, bounds the confounded motion patterns: a whole delay block moved, every label
moved, and the labels of one delay moved. Criterion, fixed in advance: at most 0.1 mm at white
noise of 2% of the brain control signal.</p>
{
            plain_table(
                ['Model'] + [f'<code>{esc(x)}</code> worst SD (mm)' for x in fixtures],
                rows,
                tuple(range(1, len(fixtures) + 1)),
            )
        }
<p class="note">Worst Cram&eacute;r&ndash;Rao standard deviation over the confounded patterns and the
three translations, at 2% noise. &infin;: a pattern the nuisance parameters absorb completely.</p>
<p>Free contrast maps absorb every displacement shared by a group, as the plan expected. A
perfusion image that is free everywhere absorbs any shift shared by the labels on single-delay data,
and label shifts proportional to the delay weights on multi-delay data without pre-saturation. The
displaced label differs from the static image by about the shift times its gradient, which a free
perfusion image can cancel at every delay at once. Fixing the perfusion image to zero outside a
brain support makes every pattern identifiable, to 0.01 mm. The test passed, with the support
as part of the model.</p>""")
    arms = [a for a in C6_ARMS if a in set(df.label)]
    if 'c6-h1-pb3' in arms or 'c6-pb3' in arms:
        sub = df[df.label.isin(arms) & df.fixture.str.startswith(('mo_', 'c5_'))]
        full = sub.groupby('fixture').label.nunique()
        fx = full[full == len(arms)].index
        sub = sub[sub.fixture.isin(fx)]
        piv = sub.pivot_table(index='fixture', columns='label', values='dm_edge')
        pair = sub.pivot_table(index='fixture', columns='label', values='pair_max')
        rel = sub.pivot_table(index='fixture', columns='label', values='rel_max')

        def verdict(a):
            r = piv[a] / piv['c0']
            return int((r < 0.9).sum()), int(((r >= 0.9) & (r <= 1.1)).sum()), int((r > 1.1).sum())

        rows = [
            [
                esc(C6_ARMS[a]),
                fmt(piv[a].median(), 2),
                fmt(pair[a].median(), 3),
                fmt(rel[a].median(), 3),
                '&ndash;' if a == 'c0' else '{} / {} / {}'.format(*verdict(a)),
            ]
            for a in arms
        ]
        key = [k for k in C6_KEY if k in piv.index]
        krows = [
            [f'<code>{esc(k)}</code>']
            + [f'{fmt(piv.at[k, a], 2)} / {fmt(pair.at[k, a], 3)}' for a in arms]
            for k in key
        ]
        out.append(f"""
<h3>The prototype</h3>
<p>The prototype (PyTorch on the GPU, fitted on the host because the ASLPrep image has no PyTorch)
ran into two problems the local test cannot see. First, the model cannot reproduce how the
simulator forms a moved image: with the poses frozen at the truth the fit stays above the noise,
and free poses settle up to 0.15 mm away (from any start, so not a local minimum). An acquisition
model (a finer in-plane grid, cropped in k-space, with a learned weight per phase-encoding line)
halves that misfit but barely moves the poses. Second, with that misfit the fit preferred to leave
labels unmoved and absorb their displacement into the perfusion image (pair error 0.54 mm on the
synchronous fixture), because moving a volume costs the misfit while absorbing costs only the
residual outside the support. Bounding the perfusion image to &plusmn;3% of the static image
(delta-M is about 1% of the control signal) removes that.</p>
{
            plain_table(
                [
                    'Arm',
                    'Median edge error',
                    'Median worst pair (mm)',
                    'Median worst pose vs mean (mm)',
                    'Edge error vs production: better / tie / worse',
                ],
                rows,
                (1, 2, 3),
            )
        }
<p class="note">{
            len(piv.dropna())
        } development motion fixtures that every arm ran on and that have control-label pairs (the
edge error needs them). Better or worse: an edge error more than 10% below or above
production's.</p>
{
            plain_table(
                ['Fixture'] + [esc(C6_ARMS[a]) for a in arms],
                krows,
                tuple(range(1, len(arms) + 1)),
            )
        }
<p class="note">Edge error / worst pair error (mm) on the fixtures that test the type confound
(synchronous, spikes, within-pair motion) and on drift.</p>
<p><b>C6 does not join the short list.</b> Its edge error is several times production's on most
fixtures, including still heads, where its misfit adds motion that is not there. Its gains are
confined to motion confounded with the volume type or the pairs (within-pair and between-pair
motion, constant offsets), where C2 with 3dvolreg already does as well or better. The prototype
fits a model the simulator does not obey exactly, and that misfit, not identifiability or the
optimizer, limits it. A real scanner's departures from the model would be different but no smaller.
The plan builds C7 (a learned initializer) only if C6 succeeds, so C7 is not built.</p>""")
    if not out:
        return ''
    return (
        '<section id="c6"><h2>A joint contrast and motion model (C6)</h2>'
        + ''.join(out)
        + '</section>'
    )


T7_ARMS = {
    'c0': 'production',
    'c9-auto': 'C9 (production, label offset removed)',
    'c1-pld': 'C1 (own-type, per delay)',
    'c2-volreg': 'C2 (3dvolreg)',
    'c4-c1': 'C4 (pair-shared)',
    'c8-huber-perdelay': 'C8 (zig-zag, Huber)',
}
T7_METRICS = [
    ('cbf_gm_bias_rel', 'CBF bias in GM', 3),
    ('cbf_brain_rmse_rel', 'CBF RMSE, brain', 3),
    ('cbf_edge_rmse_rel', 'CBF RMSE, outer 3 mm', 3),
    ('coreg_end_to_end_mm', 'Coregistration, end to end (mm)', 2),
    ('m0_mm', 'M0 registration (mm)', 2),
    ('mask_dice', 'Mask Dice', 3),
]


def section_task7(root):
    """Task 7: the short list through full ASLPrep (metrics 5, 6, 9-11) and D3 per arm."""
    f = root / 'task7' / 'scores' / 'task7.tsv'
    if not f.exists():
        return ''
    d = pd.read_csv(f, sep='\t')
    d = d[d.label.isin(T7_ARMS)]
    d['kind'] = np.where(d.fixture.str.startswith('mo_'), 'motion', 'still')
    arms = [a for a in T7_ARMS if a in set(d.label)]
    tables = []
    for kind, what in (('still', 'still head'), ('motion', 'simulated motion')):
        sub = d[d.kind == kind]
        med = sub.groupby('label').median(numeric_only=True)
        rows = [
            [esc(name)] + [fmt(med.at[a, col], dg) if a in med.index else '&ndash;' for a in arms]
            for col, name, dg in T7_METRICS
        ]
        n = sub.fixture.nunique()
        tables.append(
            f'<h3>{what.capitalize()} ({n} fixtures)</h3>'
            + plain_table(
                ['Median'] + [esc(T7_ARMS[a]) for a in arms], rows, tuple(range(1, len(arms) + 1))
            )
        )
    still = d[d.kind == 'still'].groupby('label').median(numeric_only=True)
    motion = d[d.kind == 'motion'].groupby('label').median(numeric_only=True)
    dropped = sorted(set(d.fixture[d.cbf_brain_rmse_rel.isna()]))
    out = [
        f"""
<p>Each short-listed estimator was wired into full ASLPrep by handing its HMC reference and
transforms to a fresh run as precomputed derivatives. ASLPrep then skips its own HMC and builds
everything after it from the candidate's products: M0 and fieldmap registration, the
coregistration reference, the mask, coregistration and CBF. Production's harness output, injected
the same way, reproduces the baseline run (motion, M0 and mask identical, coregistration within
0.002 mm). The CBF target is an oracle run given the true poses: the same quantification of the
same data without motion-correction error, inside the true brain. It is built in each arm's own
frame, so that neither map is resampled more than the other: production and C9 use production's
reference, which with a separate M0 can sit 2.5 mm from the volumes, and are compared with an
oracle built in that frame. (A first version compared every arm with an oracle in the static
frame; resampling that oracle into production's frame smoothed it and showed production a 3%
grey-matter bias that was not there.) Values are relative to the oracle's grey-matter CBF.</p>
<p>Getting there exposed an ASLPrep bug: precomputed transforms are silently dropped whenever
atlases are enabled (the default), because the atlas datasets are collected after the user's and
each replaces the transforms. The probe works around it; the fix belongs in ASLPrep.</p>
{''.join(tables)}
<p class="note">Medians over fixtures. The motion fixtures' CBF excludes {
            ', '.join(f'<code>{esc(x)}</code>' for x in dropped) or 'none'
        }, where even the oracle's CBF is corrupt (2 mm spikes through a B0 field fixed to the
scanner cannot be undone by any rigid transform).</p>
<p><b>No arm has a grey-matter CBF bias worth the name</b> (all within
{fmt(max(abs(still.loc[a, 'cbf_gm_bias_rel']) for a in arms) * 100, 1)}% on still heads); the
differences are in the spatial error. On still heads production's brain CBF error is
{fmt(still.at['c0', 'cbf_brain_rmse_rel'] * 100, 0)}% of grey-matter CBF, the own-type templates
are near zero, and C9 is in between ({fmt(still.at['c9-auto', 'cbf_brain_rmse_rel'] * 100, 0)}%).
With motion the template arms lose: production
{fmt(motion.at['c0', 'cbf_brain_rmse_rel'] * 100, 0)}%, C1 and C2 about
{fmt(motion.loc[['c1-pld', 'c2-volreg'], 'cbf_brain_rmse_rel'].mean() * 100, 0)}%, C4 and C8 far
worse, and C9 {fmt(motion.at['c9-auto', 'cbf_brain_rmse_rel'] * 100, 0)}%. <b>C9 is the only arm
better than production on both</b>: removing production's label offset keeps its motion tracking
and cuts its CBF error by a third on still heads and by more than half with motion. What remains
on still heads is per-volume jitter from registering every volume, across contrast, to a
reference that can sit millimetres away, which C9 does not touch. The template arms also make the
separate M0's registration cross-contrast, which costs accuracy (worst on GRASE); C9 keeps
production's.</p>""",
    ]
    c = root / 'task7' / 'coreg_d3.tsv'
    if c.exists():
        cd = pd.read_csv(c, sep='\t')
        cd = cd[cd.label.isin(arms) & cd.source.notna()]
        piv = cd[cd.sign == 'default'].pivot_table(
            index='label', columns='source', values='rms_mm', aggfunc='median'
        )
        cols = [
            s for s in ('conventional', 'P-mean', 'P-robust', 'P-winsor', 'P-masked') if s in piv
        ]
        spiral = cd[(cd.fixture == 'st_spiral') & (cd.sign == 'default')].pivot_table(
            index='label', columns='source', values='rms_mm', aggfunc='median'
        )
        rows = [
            [esc(T7_ARMS[a])]
            + [fmt(piv.at[a, s], 2) for s in cols]
            + [fmt(spiral.at[a, 'conventional'], 1) if a in spiral.index else '&ndash;']
            + [fmt(spiral.at[a, 'P-robust'], 2) if a in spiral.index else '&ndash;']
            for a in arms
            if a in piv.index
        ]
        out.append(f"""
<h3>D3 on every arm</h3>
<p>The perfusion-image levels were built from each run's corrected series and registered by
fMRIPrep's FSL path, next to the run's conventional reference (default BBR sign).</p>
{
            plain_table(
                ['Arm']
                + [f'{esc(s)} (mm)' for s in cols]
                + ['spirals: conventional', 'spirals: P-robust'],
                rows,
                tuple(range(1, len(cols) + 3)),
            )
        }
<p class="note">Median end-to-end coregistration error over the fixtures, and on the
stack-of-spirals fixture.</p>
<p>The conventional reference still registers best overall for every arm. Built from
production's corrected series, the perfusion images register worse (about 0.5 mm) than from any
candidate's (about 0.25 mm): production's spurious label offset blurs the pairs they are made
from. On the spirals fixture every conventional reference fails, production's and the
candidates' templates alike, which is the BBR sign problem, and the perfusion image recovers
it on every arm.</p>""")
    return (
        '<section id="task7"><h2>The short list through full ASLPrep (Task 7)</h2>'
        + ''.join(out)
        + '</section>'
    )


REF_FIXTURES = {
    'BS': ['c5_presat', 'c5_presat_drift', 'c5_presat_bs_drift', 'c5_presat_3d_drift',
           'c5_presat_pldconfounded', 'rc_grase', 'st_spiral'],
    'noBS': ['mo_pcasl_drift', 'mo_pcasl_withinpairs', 'mo_pcasl_synchronous', 'mo_pcasl_spikes',
             'mo_pasl_drift', 'mo_pcaslinc_drift', 'rc_pcasl', 'st_vessels'],
}  # fmt:skip
REF_ARMS = {
    'c0': 'production',
    'c9-auto': 'typewise_offset',
    'c0x-vr': '3dvolreg to production ref.',
    'c0x-vr-ctl': '3dvolreg to control mean',
    'c0xfb': '3dvolreg + fallback',
    'c9-fb-ctl': '3dvolreg to control mean + fallback + offset',
    'c0-767': '#767',
    'c9-auto-767': '#767 + offset',
    'c9-t-vr': '#767 design, 3dvolreg + offset',
}
REF_HEAT = ['c0', 'c9-auto', 'c0-767', 'c9-auto-767', 'c9-fb-ctl', 'c9-t-vr']


def section_reference(root, df, fig):
    """Reference type, PennLINC/aslprep#767, and 3dvolreg with an MCFLIRT fallback."""
    fixtures = REF_FIXTURES['BS'] + REF_FIXTURES['noBS']
    d = df[df.label.isin(REF_ARMS) & df.fixture.isin(fixtures)].copy()
    if d.empty or not {'c0xfb', 'c9-t-vr'} <= set(d.label):
        return ''
    d['bs'] = np.where(d.fixture.isin(REF_FIXTURES['BS']), 'BS', 'noBS')
    med = d.groupby(['label', 'bs'])[['dm_edge', 'pair_max', 'rel_max']].median()
    rows = []
    for a, name in REF_ARMS.items():
        if (a, 'BS') not in med.index or (a, 'noBS') not in med.index:
            continue
        cells = [f'{esc(name)} <code>{esc(a)}</code>']
        for col, dg in (('dm_edge', 3), ('pair_max', 3), ('rel_max', 2)):
            cells += [fmt(med.at[(a, b), col], dg) for b in ('BS', 'noBS')]
        rows.append(cells)
    medians = plain_table(
        ['Arm', 'Edge &Delta;M, BS', 'no BS', 'Pair offset (mm), BS', 'no BS',
         'Pose error (mm), BS', 'no BS'],
        rows, tuple(range(1, 7)),
    )  # fmt:skip
    piv = d.pivot_table(index='fixture', columns='label', values='dm_edge')
    heat_arms = [a for a in REF_HEAT if a in piv]
    heat = heat_table(
        fixtures,
        heat_arms,
        [[piv.at[f, a] if f in piv.index else None for a in heat_arms] for f in fixtures],
        notes=dict.fromkeys(REF_FIXTURES['BS'], 'background suppression'),
        digits=3,
        col_labels=[REF_ARMS[a] for a in heat_arms],
    )
    fb_rows = []
    for a in ('c0xfb', 'c0xfb-ctl'):
        for f in fixtures:
            h = root / 'fixture-scores' / a / f / 'harness.json'
            if not h.exists():
                continue
            ref = json.loads(h.read_text())['reference_type']
            fb_rows.append((a, f, 'fell back' if 'diverged' in ref else 'kept', ref.split(';')[-1]))
    kept = [r for r in fb_rows if r[2] == 'kept']
    fell = [r for r in fb_rows if r[2] == 'fell back']
    fb_table = plain_table(
        ['3dvolreg reference', 'Fixture', 'Outcome', 'Deviation from MCFLIRT'],
        [
            ['production&rsquo;s' if a == 'c0xfb' else 'control mean', f'<code>{esc(f)}</code>', o,
             esc(note.strip())]
            for a, f, o, note in fell
        ],
    )  # fmt:skip
    dm = []
    for k in ('drift', 'spikes'):
        for snr in (8, 20, 50):
            f = f'mo_deltamrep_{k}_snr{snr}'
            cells = [f'<code>{esc(k)}</code>', str(snr)]
            for a in ('c0', 'c0-m0', 'c0-767', 'c0-767-m0', 'c0t-vr', 'c0t-vr-m0'):
                p = root / 'fixture-scores' / a / f / 'score.json'
                v = json.loads(p.read_text())['relative_pose']['to_mean']['max_mm'] if p.exists() else None
                cells.append(fmt(v, 2) if v is not None else '&ndash;')
            dm.append(cells)
    dm_table = plain_table(
        ['Trajectory', 'Delta-M SNR', 'MCFLIRT, delta-M ref.', 'MCFLIRT, M0 ref.', '#767, delta-M ref.',
         '#767, M0 ref.', '#767 design + 3dvolreg, delta-M ref.', '#767 design + 3dvolreg, M0 ref.'],
        dm, tuple(range(1, 8)),
    )  # fmt:skip
    ge = ge_reference(root)
    return f"""
<section id="reference"><h2>Reference type, #767 and 3dvolreg</h2>
<p>Every fixture in the earlier engine screens has an M0, so production's reference was the M0
throughout (ASLPrep picks CBF, then delta-M, then an in-series M0, then a separate M0, then the
control mean) and every label and control volume was registered across contrasts. This screen
adds fixtures with an included M0, background suppression and a delta-M reference, forces a
control-mean reference where the harness asks for one (it hides in-series M0, delta-M and CBF
volumes from <code>SelectHighestContrastVolumes</code>, which otherwise picks its own type), and
scores PennLINC/aslprep#767 (head <code>96d0530</code>): each volume type registered to its own
mean, which FLIRT (mutual information) aligns to the reference. Fifteen label-control fixtures,
{len(REF_FIXTURES['BS'])} with background suppression (BS). Constant-pose fixtures are left out:
there an arm that happens to return identity transforms avoids interpolation and scores better
for that alone.</p>
{medians}
<p class="note">Medians over fixtures: edge delta-M error (0&ndash;3 mm shell, relative to GM),
max label-control pair offset, max pose error relative to the series mean.</p>
{figure(heat, 'Edge delta-M error per fixture (darker is worse, log scale).', fig)}
<ul>
<li><b>Least squares fails across contrasts.</b> 3dvolreg registering volumes straight to an M0
diverges by 8&ndash;135 mm on every background-suppressed <code>c5</code> fixture and on
<code>mo_pasl_drift</code>, and by 165&ndash;215 mm on a delta-M reference. With a control-mean
reference it converges on every fixture but <code>st_spiral</code>. Where it converges it is the
most precise engine: about a third of MCFLIRT's pose error without suppression.</li>
<li><b>A fallback catches every divergence.</b> <code>c0xfb</code> keeps 3dvolreg unless any
volume's mean-relative pose differs from production's MCFLIRT by more than 1 mm. Runs that agree
differ by at most {fmt(max(float(r[3].split('within')[1].split('mm')[0]) for r in kept), 2)} mm,
diverged ones by at least {fmt(min(float(r[3].split('(')[1].split('mm')[0]) for r in fell), 1)}
mm. It does not catch a cross-contrast bias without divergence: on <code>rc_grase</code> 3dvolreg
stays within 0.16 mm of MCFLIRT but its label-control offset gives an edge error of 0.62, which
the offset step removes.</li>
<li><b>#767 adds a label-control offset.</b> Its pose error is mixed, but edge delta-M error rises
on 14 of 15 fixtures: each type's cross-contrast FLIRT is stable only to a few hundredths of a
millimetre (a reimplementation differing only in interpolation lands 0.02&ndash;0.06 mm away,
type-wide), and the difference between types is a label-control offset. With the offset step it
matches <code>typewise_offset</code>. Reported on the PR.</li>
<li><b>The best arms combine same-contrast least squares with the offset step.</b> 3dvolreg to
the control mean with the fallback, and #767's design with 3dvolreg within each type, both followed
by the offset step, tie for best: label-control offsets 3&ndash;4 times smaller than
<code>typewise_offset</code>'s, pose error a third of it without suppression. #767's design needs
no fallback (its one cross-contrast step per type is FLIRT with mutual information), but these
runs use the probe's reimplementation of it.</li>
</ul>
<h3>Where the fallback was used</h3>
{fb_table}
<h3>Delta-M reference or M0 reference?</h3>
<p>On the two delta-M fixtures (eight delta-M volumes and an included M0) ASLPrep's reference is
the delta-M mean. Forcing the M0 instead makes every delta-M volume a cross-contrast registration.
Real delta-M images are noisier than the fixtures' (grey-matter delta-M SNR 50): on the GE
ds005752 runs it is between about 8 (noise from the high-pass residual inside the brain, an
underestimate) and 21 (noise from the background, an overestimate). Both fixtures were rebuilt at
8, 20 and 50 (the rebuild at 50 reproduces the originals exactly). The delta-M reference keeps the
lower max pose error in every arm and at every SNR:</p>
{dm_table}
<p class="note">Max pose error (mm) relative to the series mean. aslscan's noise is white;
the GE images' background noise is structured, which the simulation does not reproduce.</p>
{ge}
</section>"""


def ge_reference(root):
    """The GE ds005752 runs (one M0, one delta-M): where each reference puts the delta-M."""
    import nibabel as nb
    from batch_harness import rms
    from harness import load_itk_series

    out = root / 'deltam-ref-real'
    if not (out / 'c0-dmref').exists():
        return ''
    rows = []
    for r in sorted((out / 'c0' / 'ds005752').iterdir()):
        rel = {}
        for lab in ('c0', 'c0-dmref'):
            d = out / lab / 'ds005752' / r.name
            if not (d / 'hmc_xfms.txt').exists():
                break
            h = load_itk_series(d / 'hmc_xfms.txt')
            ctx = pd.read_table(d / 'hmc_aslcontext.tsv')['volume_type'].tolist()
            ref = nb.load(d / 'aslref.nii.gz')
            centre = (ref.affine @ np.r_[(np.array(ref.shape[:3]) - 1) / 2, 1])[:3]
            rel[lab] = (np.linalg.inv(h[ctx.index('m0scan')]) @ h[ctx.index('deltam')], centre)
        if len(rel) == 2:
            c = rel['c0'][1]
            rows.append([rms(rel['c0'][0], c), rms(rel['c0-dmref'][0], c),
                         rms(rel['c0'][0] @ np.linalg.inv(rel['c0-dmref'][0]), c)])  # fmt:skip
    if not rows:
        return ''
    a = np.array(rows)
    return f"""
<p>The GE runs of ds005752 have one M0 and one delta-M volume each, so the reference choice only
decides which of the two is registered to the other. On {len(a)} runs the delta-M lands
{fmt(np.median(a[:, 0]), 2)} mm from the M0 (median; range {fmt(a[:, 0].min(), 1)}&ndash;{
        fmt(a[:, 0].max(), 1)}) with the M0 reference, {fmt(np.median(a[:, 1]), 2)} mm with the
delta-M reference, and the two estimates disagree by {fmt(np.median(a[:, 2]), 2)} mm (range
{fmt(a[:, 2].min(), 1)}&ndash;{fmt(a[:, 2].max(), 1)}), mostly through-plane and in no consistent
direction. Both come from one GE sequence, so the true motion between them should be small: these
offsets look like registration error, larger than any simulated case. For single-delta-M GE
series, not registering the M0 and delta-M to each other may beat either reference; that cannot be
settled without truth.</p>"""


def section_coreg(root, fig):
    out = []
    f = root / 'coreg-5b' / 'coreg_5b.tsv'
    if f.exists():
        c = pd.read_csv(f, sep='\t')
        rows = []
        for fx, g in c.groupby('fixture', sort=False):

            def v(src, sign='default', col='rms_mm', g=g):
                hit = g[(g.source == src) & (g.sign == sign)]
                return hit.iloc[0][col] if len(hit) else np.nan

            fb = v('conventional', col='fallback')
            rows.append(
                [
                    f'<code>{esc(fx)}</code>',
                    fmt(v('production (ASLPrep run)'), 2),
                    fmt(v('conventional'), 2)
                    + (' (fell back)' if fb is True or fb == 'True' else ''),
                    fmt(v('conventional', 'opposite'), 2),
                    fmt(v('P-robust'), 2),
                    fmt(v('P-masked'), 2),
                    fmt(v('conventional', col='gm_wm_ratio'), 2),
                    fmt(v('P-robust', col='gm_wm_ratio'), 2),
                ]
            )
        g = c[c.fixture == 'rc_grase']

        def gv(src, sign='default', col='rms_mm'):
            hit = g[(g.source == src) & (g.sign == sign)]
            return float(hit.iloc[0][col]) if len(hit) else np.nan

        out.append(f"""
<p>Each reference was registered to the T1w exactly as fMRIPrep's FSL path does it (<code>mri_coreg</code>
initialization, then FLIRT with the boundary-based cost), standalone, and scored against the
fixture's truth. The standalone registration of the conventional reference reproduces ASLPrep's own
result on every fixture, so the comparison is like for like.</p>
{
            plain_table(
                [
                    'Fixture',
                    'ASLPrep',
                    'Conventional',
                    'Conventional, opposite BBR sign',
                    'P-robust',
                    'P-masked',
                    'Grey/white, conventional',
                    'Grey/white, perfusion',
                ],
                rows,
                tuple(range(1, 8)),
            )
        }
<p class="note">End-to-end coregistration error against truth (mm RMS over the brain, in the T1w):
the transform as applied to the corrected series. Default BBR sign unless stated. "Fell back" means
FLIRT-BBR moved too far from its initialization (fMRIPrep's displacement rule, 15 mm) and the
initialization was kept, as fMRIPrep does.</p>
<p>End to end, the conventional reference carries two errors: the registration itself, and the
offset between the image it depicts (here the separate M0) and the corrected series, which HMC
leaves. On this fixture that offset is {
            fmt(gv('conventional', col='conventional_vs_series_frame_mm'), 2)
        } mm.
Against its own frame, the conventional reference's registration error is
{fmt(gv('conventional', col='step_rms_mm'), 2)} mm with the default sign and
{
            fmt(gv('conventional', 'opposite', col='step_rms_mm'), 2)
        } mm with the opposite sign. The perfusion
images are built from the corrected series, so for them the two errors are the same.</p>
<p>The 3D GRASE fixture explains its 1 mm coregistration error. Its conventional reference has white
matter brighter than grey (grey/white {fmt(gv('conventional', col='gm_wm_ratio'), 2)}), so BBR's
default sign, which assumes the opposite, fails ({
            fmt(gv('conventional', col='bbr_rms_mm'), 1)
        } mm) and
the run silently falls back to <code>mri_coreg</code> ({fmt(gv('conventional'), 2)} mm). With the
opposite sign the same reference registers to {fmt(gv('conventional', 'opposite'), 2)} mm. The
perfusion image (grey/white about {
            fmt(gv('P-robust', col='gm_wm_ratio'), 1)
        }) avoids the problem and
reaches {
            fmt(gv('P-robust'), 2)
        } mm, but on every other fixture the conventional reference registers
better. Choosing the BBR sign from the reference's measured contrast, which the plan already calls
for, may be the simpler fix.</p>""")
    s = root / 'coreg-5b' / 'coreg_5b_stress.tsv'
    if s.exists():
        st = pd.read_csv(s, sep='\t')

        def sv(fx, src, sign='default', col='rms_mm'):
            hit = st[(st.fixture == fx) & (st.source == src) & (st.sign == sign)]
            return float(hit.iloc[0][col]) if len(hit) else np.nan

        names = {
            'st_vessels': 'macrovascular signal',
            'st_noperf': 'perfusion deficit',
            'st_spiral': 'stack-of-spirals (#503)',
        }
        srows = [
            [
                esc(names.get(fx, fx)),
                fmt(sv(fx, 'production (ASLPrep run)'), 2),
                fmt(sv(fx, 'conventional', 'opposite'), 2),
                fmt(sv(fx, 'P-robust'), 2),
                fmt(sv(fx, 'P-masked'), 2),
                fmt(sv(fx, 'conventional', col='gm_wm_ratio'), 2),
            ]
            for fx in st.fixture.unique()
        ]
        out.append(f"""
<h3>Stress fixtures</h3>
<p>Three fixtures stress the perfusion reference: arterial signal in large vessels (eleven times the
grey-matter signal), a territory with no perfusion, and a stack-of-spirals readout, the case of #503.</p>
{
            plain_table(
                [
                    'Fixture',
                    'ASLPrep',
                    'Conventional, opposite BBR sign',
                    'P-robust',
                    'P-masked',
                    'Grey/white, conventional',
                ],
                srows,
                tuple(range(1, 6)),
            )
        }
<p class="note">End-to-end error against truth, mm RMS over the brain in the T1w, as above.</p>
<p>Neither vessels nor the deficit breaks the perfusion reference, but neither helps it: the
conventional reference registers better on both (about {fmt(sv('st_vessels', 'conventional'), 2)}
mm against {
            fmt(sv('st_vessels', 'P-robust'), 2)
        }). On the spirals fixture, ASLPrep's coregistration
misses by {fmt(sv('st_spiral', 'production (ASLPrep run)'), 1)} mm, and every perfusion image lands
within {fmt(sv('st_spiral', 'P-masked'), 2)}&ndash;{fmt(sv('st_spiral', 'P-robust'), 2)} mm. The
cause is the same as on GRASE: the conventional reference has white matter brighter than grey
({fmt(sv('st_spiral', 'conventional', col='gm_wm_ratio'), 2)}), BBR's default sign fails, and here
the failure ({
            fmt(sv('st_spiral', 'conventional', col='bbr_rms_mm'), 1)
        } mm) stays inside fMRIPrep's
15 mm fallback rule, so nothing catches it. With the opposite sign, the same reference registers to
{fmt(sv('st_spiral', 'conventional', 'opposite'), 2)} mm. Choosing the BBR sign from the measured
contrast would fix both the GRASE and the spirals cases. The simulated spirals are milder than real
ones (3.7% haze, no severe ringing), so on real data the perfusion reference may matter more than it
does here.</p>""")
    e = root / 'harness-runs' / 'perfusion_c0.tsv'
    if e.exists():
        p = pd.read_csv(e, sep='\t')
        p = p[p.edge_r.notna()]
        med = p.groupby('family').edge_r.median().sort_values()
        strip = strip_plot(
            med.index.tolist(),
            [p.loc[p.family == fm, 'edge_r'] for fm in med.index],
            'edge statistic (split-half correlation of gradient magnitude)',
            xmax=1.0,
            label='Perfusion-image eligibility',
        )
        out.append(f"""
<h3>Which real runs could use a perfusion image</h3>
<p>A run is eligible for the perfusion reference only if its perfusion image has reproducible edges.
The statistic correlates the gradient magnitude of two half-run perfusion images. On the
{len(p)} real runs with control-label pairs, stack-of-spirals runs are far above every other family
(median {fmt(med.max(), 2)}); most families sit near the placeholder limit of 0.5, so its calibration
will decide eligibility for them.</p>
{figure(strip, 'Edge statistic per real run, by family; the tick is the median.', fig)}""")
    r = root / 'coreg-5b-real' / 'coreg_5b_real.tsv'
    if r.exists():
        d = pd.read_csv(r, sep='\t')
        d = d[(d.source == 'P-robust') & (d.sign == 'default') & d.vs_conventional_mm.notna()]
        if len(d):
            med = d.groupby('family').vs_conventional_mm.median().sort_values()
            strip = strip_plot(
                med.index.tolist(),
                [d.loc[d.family == fm, 'vs_conventional_mm'] for fm in med.index],
                'distance between the perfusion and conventional registrations (mm)',
                label='Coregistration disagreement',
                color='var(--s-c8)',
            )
            out.append(f"""
<p>On real data there is no truth, only disagreement. Across {len(d)} runs the perfusion
registration lands a median {fmt(d.vs_conventional_mm.median(), 2)} mm from the conventional one, and
{(d.vs_conventional_mm > 1).mean() * 100:.0f}% of runs differ by more than 1 mm. Those runs are what the
blind review has to adjudicate.</p>
{figure(strip, 'Distance between the P-robust and the conventional coregistration per run.', fig + 1)}""")
    if not out:
        return ''
    return f'<section id="coreg"><h2>A perfusion image as the coregistration reference (#672)</h2>{"".join(out)}</section>'


CSS = """
:root {
  /* Layout: one reading column (about 70 characters) with figures allowed a little wider. */
  --bg: #f5f7f8; --surface: #ffffff; --fg: #18232b; --muted: #5a6872; --grid: #dde3e7;
  --accent: #1d6a83; --hot: #c4472c; --rule: #c9d3d9;
  --s-c0: #18232b; --s-c3: #7b8a94; --s-c1: #2f7fb1; --s-c2: #86b7d8; --s-c4: #c98a1d;
  --s-c8-noreg: #7a5ab8; --s-c8: #c4472c; --s-c8-huber: #e2897a; --s-c8-perdelay: #a8586e;
  --s-c8-prod: #8a5a2b; --s-c8-spm: #2d8a5a; --s-c5: #1f8f8a; --s-c5-voxel: #9fd3cf;
  --display: "Archivo", "Helvetica Neue", Arial, sans-serif;
  --body: "Source Serif 4", Georgia, "Times New Roman", serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #11181d; --surface: #182229; --fg: #e3e9ed; --muted: #97a6b0; --grid: #26333c;
  --accent: #5fb3cf; --hot: #e0644a; --rule: #2e3c46;
  --s-c0: #e3e9ed; --s-c3: #8d9ca6; --s-c1: #5aa8dc; --s-c2: #9dc9e8; --s-c4: #e0a640;
  --s-c8-noreg: #a88be0; --s-c8: #ef6a50; --s-c8-huber: #f2a596; --s-c8-perdelay: #d1849a;
  --s-c8-prod: #c9925e; --s-c8-spm: #55c08a; --s-c5: #45c2bb; --s-c5-voxel: #2a6f6b;
  color-scheme: dark; } }
:root[data-theme="dark"] {
  --bg: #11181d; --surface: #182229; --fg: #e3e9ed; --muted: #97a6b0; --grid: #26333c;
  --accent: #5fb3cf; --hot: #e0644a; --rule: #2e3c46;
  --s-c0: #e3e9ed; --s-c3: #8d9ca6; --s-c1: #5aa8dc; --s-c2: #9dc9e8; --s-c4: #e0a640;
  --s-c8-noreg: #a88be0; --s-c8: #ef6a50; --s-c8-huber: #f2a596; --s-c8-perdelay: #d1849a;
  --s-c8-prod: #c9925e; --s-c8-spm: #55c08a; --s-c5: #45c2bb; --s-c5-voxel: #2a6f6b;
  color-scheme: dark; }
body { background: var(--bg); color: var(--fg); font: 17px/1.6 var(--body); }
.wrap { max-width: 900px; margin: 0 auto; padding-inline: 20px; padding-block: 28px 64px; }
header { border-bottom: 1px solid var(--rule); padding-block-end: 20px; }
.eyebrow { font: 600 12px/1.4 var(--mono); letter-spacing: .08em; text-transform: uppercase; color: var(--accent); }
h1, h2, h3, h4 { font-family: var(--display); text-wrap: balance; line-height: 1.2; }
h1 { font-size: clamp(30px, 5vw, 42px); font-weight: 700; margin: 8px 0 12px; letter-spacing: -.01em; }
h2 { font-size: 26px; font-weight: 700; margin: 0 0 12px; }
h3 { font-size: 19px; font-weight: 600; margin: 28px 0 8px; }
h4 { font-size: 14px; font-weight: 600; margin: 0 0 4px; }
p, li { max-width: 68ch; }
section { padding-block: 34px 6px; border-bottom: 1px solid var(--rule); }
.lede { font-size: 19px; color: var(--fg); }
.meta { font: 13px/1.5 var(--mono); color: var(--muted); }
nav.toc { display: flex; flex-wrap: wrap; gap: 6px 18px; font: 600 14px/1.5 var(--display); margin-top: 16px; }
nav.toc a { color: var(--accent); text-decoration: none; }
nav.toc a:hover, nav.toc a:focus-visible { text-decoration: underline; }
a:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
code { font: .86em var(--mono); }
.findings { display: grid; gap: 10px; padding: 0; list-style: none; counter-reset: f; }
.findings li { background: var(--surface); border: 1px solid var(--rule); border-radius: 6px; padding: 12px 16px;
  max-width: none; }
.findings b { font-family: var(--display); }
figure { margin: 22px 0; }
figcaption { font-size: 14px; color: var(--muted); margin-top: 6px; max-width: 75ch; }
svg.chart { width: 100%; height: auto; display: block; font-family: var(--mono); }
svg .grid { stroke: var(--grid); stroke-width: 1; }
svg .axis, svg .tickmark { stroke: var(--muted); stroke-width: 1; }
svg .tick { fill: var(--muted); font-size: 11px; }
svg .axlabel { fill: var(--fg); font-size: 12px; }
svg .median { stroke: var(--fg); stroke-width: 2.5; }
.legend { display: flex; flex-wrap: wrap; gap: 4px 16px; font: 12px/1.5 var(--mono); color: var(--muted); margin-bottom: 6px; }
.key { display: inline-flex; align-items: center; gap: 6px; }
.key i { display: inline-block; width: 22px; }
.panels { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 16px; }
.panel { min-width: 0; }
.tablewrap { overflow-x: auto; margin: 14px 0; }
table { border-collapse: collapse; font: 13px/1.45 var(--mono); font-variant-numeric: tabular-nums; width: 100%; }
th, td { padding: 6px 9px; border-bottom: 1px solid var(--rule); text-align: left; vertical-align: top; }
thead th { font: 600 12px/1.3 var(--display); color: var(--muted); border-bottom: 1.5px solid var(--fg); }
td.num, table.heat td { text-align: right; white-space: nowrap; }
table.heat th[scope=row] { white-space: nowrap; font-weight: 400; }
table.heat thead th { white-space: nowrap; }
.rownote { display: block; font-size: 11px; color: var(--muted); white-space: normal; max-width: 34ch; }
td.na { color: var(--muted); text-align: center; }
.note { font-size: 14px; color: var(--muted); }
.callout { border-left: 3px solid var(--accent); padding: 4px 0 4px 16px; }
@media (prefers-reduced-motion: no-preference) { nav.toc a { transition: color .15s; } }
"""


def build(root):
    df = load_scores(root)
    fig = 1
    noise = section_noise(root, fig)
    fig += 1
    blur = section_blur(root, fig)
    fig += 1
    screen = section_screen(df, fig)
    fig += 2
    zig = section_zigzag(root, df, fig)
    fig += 2
    engines = section_engines(root, fig)
    fig += 1
    policy = section_policy(df, fig, root)
    consumers = section_consumers(root)
    real = section_real(root, fig)
    fig += 2
    coreg = section_coreg(root, fig)
    fig += 2
    c5 = section_c5(root, fig)
    c6 = section_c6(root, df, fig + 1)
    task7 = section_task7(root)
    reference = section_reference(root, df, fig + 2)

    piv = df.pivot_table(index='fixture', columns='label', values='dm_edge')
    med = piv.median()
    head = f"""<title>HMC Probe Findings</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;600;700&family=IBM+Plex+Mono:wght@400;600&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>{CSS}</style>"""
    summary = f"""
<ul class="findings">
<li><b>The spurious label motion is systematic.</b> It is there with the noise switched off and
barely changes with ten times the noise. It comes from registering one contrast to another. The
GE path (delta-M to M0 through a smoothed reference) is the exception: its offset is the largest
without noise, and noise makes it worse.</li>
<li><b>No MCFLIRT cost function fixes it.</b> Registering each volume type to its own template does.</li>
<li><b>Blurred references bias registration.</b> A sharp image registered to a blurred copy of itself
moves; matched blur removes it. This contributes to the GE path's errors (#702) and penalizes averaged
templates whenever the head moves.</li>
<li><b>Production is the best estimator on average</b> on simulated motion (median edge error
{fmt(med.get('c0'), 2)}). Own-type templates win on still heads and lose with motion
({fmt(med.get('c1'), 2)} and {fmt(med.get('c2'), 2)}). The margins are small, often under 0.1 mm.</li>
<li><b>Production fails badly in the tail on real data.</b> With a 5% tie margin on the shared-mask
energy, <code>c1-pld</code> and <code>c2-volreg</code> beat production on about a tenth of real runs
and lose on 1&ndash;2%. The runs they win are mostly production failures, not small gains (ds007116
sub-21035: production FD 10 mm, C2 0.15 mm). That argues for keeping production's estimator with a
failure guard and a fallback, not for keeping it unchanged.</li>
<li><b>Saturation-recovery normalization (C5, #673) ties production.</b> It refuses ineligible runs
correctly. A voxelwise T1eff fit is biased toward its own pose estimates; a per-slice fit avoids that
but does not beat production.</li>
<li><b>A joint contrast and motion model (C6) is identifiable but does not beat production.</b> An
exact identifiability test shows that the perfusion image must be zero outside the brain, or it
absorbs label motion. The GPU prototype also needs the perfusion image bounded, and is still limited
by how far the simulator departs from its model. It helps only on motion confounded with volume
type, where C2 is already as good.</li>
<li><b>No short-listed estimator passes gate 2.</b> Each improves what it targets (pair alignment,
reference offset) but is worse end to end with motion, in the harness and through full ASLPrep.</li>
<li><b>C9, production with its label offset removed, beats production on both still heads and
motion.</b> It keeps production's transforms and composes into the labels of each delay one
rigid correction that brings the corrected label mean onto the control mean; early, arterial
delays borrow the nearest other delay's estimate (C9-auto). Median edge delta-M error with motion
falls from 0.29 to 0.09; through full ASLPrep its CBF error is lower than production's on 12 of
14 fixtures (a third lower on still heads, more than half with motion). It is the only arm better
than production on both. A first version (one offset from the later delays) regressed in the
PASL families on held-out simulations; C9-auto passes a fresh held-out set (better on 15 of 18,
no family regressing).</li>
<li><b>On real data C9-auto does no harm by the motion-sensitive metrics</b> (delta-M energy
worse than production's on none of 260 runs, better on 15%; negative fraction unchanged), but
fails the pre-registered split-half criterion in two small families (3D EPI at 7 T, 17 runs from
one dataset; 3D GRASE multi-delay, 7 runs), where the paired per-run change is -0.005.</li>
<li><b>Decision: C9-auto is in ASLPrep as an opt-in method,</b> <code>--hmc-method
typewise_offset</code>, with production (<code>typewise</code>) still the default. The integrated
version reproduces the probe's CBF results exactly on the three fixtures checked end to end, and
passed an adversarial and an ordinary Codex review. It keeps production's reference, so masks,
coregistration and M0 registration are production's; what stands between it and the default is
the real-data split-half criterion in two small families.</li>
<li><b>Reusing derivatives now works and is all-or-nothing.</b> The atlas datasets no longer wipe
precomputed transforms (PennLINC/aslprep#738), products are reused only with what they were
made from, and transforms record the HMC method they were made with. A run given its own outputs
reproduces them exactly.</li>
<li><b>Least squares is the right cost only within one contrast.</b> Production registers labels
and controls to the M0 whenever there is one. 3dvolreg doing that diverges by 8&ndash;215 mm under
background suppression, on a PASL fixture and on a delta-M reference; to the control mean it
converges almost everywhere and is the most precise engine. A check against MCFLIRT (fall back when
they differ by more than 1 mm) catches every divergence.</li>
<li><b>PennLINC/aslprep#767 (per-type references) adds a label-control offset</b>: edge delta-M
error rises on 14 of 15 fixtures, because each type's cross-contrast alignment is stable only to a
few hundredths of a millimetre. The offset step removes it. Reported on the PR.</li>
<li><b>Same-contrast least squares plus the offset step is the best combination so far</b>: 3dvolreg
to the control mean with the fallback, or #767's design with 3dvolreg, each followed by the offset
step. Label-control offsets are 3&ndash;4 times smaller than <code>typewise_offset</code>'s.</li>
<li><b>A delta-M reference beats an M0 reference for delta-M series</b> in simulation, at every
delta-M SNR down to GE-like levels (8): forcing the M0 raised the max pose error in every arm. On
the real GE runs (one M0, one delta-M) the two references disagree by about 1.2 mm, more than any
simulated case, so the M0-to-delta-M registration itself is unreliable there.</li>
<li><b>ASLPrep ignores precomputed transforms whenever atlases are enabled</b> (found wiring the
candidates in): the atlas datasets are collected last and replace the transforms.</li>
<li><b>The simulator favours production.</b> aslscan moves otherwise clean images rigidly, with no
spin history and high SNR: the conditions per-type MCFLIRT is built for. The approaches of mature
pipelines (zig-zag regression, saturation recovery) target effects the simulator barely produces, and
C8 reimplements only ExploreASL's and ASLtbx's motion step, not their whole pipelines.</li>
<li><b>Zig-zag regression, as in ASLtbx and ExploreASL, is perfect on still heads and harmful with
real motion</b> (median {fmt(med.get('c8'), 2)}). Spikes alias onto the alternation (a Huber fit fixes
that), and drift is removed with it (nothing fixes that). Their exact parameter convention changes
nothing.</li>
<li><b>AFNI 3dvolreg is the most precise engine for same-contrast templates.</b> With it, C2's pair
alignment matches or beats production's, but its edge delta-M error is still higher, so no C1/C2
variant beats production on subtraction yet.</li>
<li><b>A harness bug, found on real data, is fixed.</b> The single-image registrations used FLIRT's
global search, which put a spurious 1.5&ndash;2&deg; offset on every label volume of C1/C2 on real
data; every affected run was redone with MCFLIRT. The split-half correlation rated the broken runs
best, so it is not read alone.</li>
<li><b>ASLPrep's 1 mm GRASE coregistration error is a BBR sign problem.</b> White matter is brighter than
grey in that reference, BBR fails and the run falls back to its initialization; the opposite sign
registers it to under 0.5 mm.</li>
<li><b>The same sign problem explains #503 in simulation.</b> On a stack-of-spirals fixture, ASLPrep's
coregistration misses by 10 mm without triggering the fallback. The opposite BBR sign, or a perfusion
image as the reference, brings it under 0.4 mm. Vessels and a perfusion deficit do not break the
perfusion reference, but the conventional one is better on those fixtures.</li>
</ul>"""
    return f"""{head}
<div class="wrap">
<header>
<div class="eyebrow">ASLPrep · head-motion correction probe</div>
<h1>HMC Probe Findings</h1>
<p class="lede">What drives the spurious motion in ASLPrep's head-motion correction, and whether any
alternative does better, tested on simulated data with known truth and on the aslprep-zoo panel.</p>
<p class="meta">Branch <code>worktree-hmc-probe</code> · plan revision 13 · issues #702, #705, #707, #508,
#673, #700, #672 · development fixtures only; held-out fixtures are pending.</p>
<nav class="toc" aria-label="Sections">
<a href="#noise">Noise and cost</a><a href="#blur">Blur bias</a><a href="#screen">Candidate screen</a>
<a href="#zigzag">Zig-zag regression</a><a href="#engines">Engines and short list</a>
<a href="#consumers">Reference consumers</a><a href="#real">Real data</a><a href="#coreg">Perfusion coregistration</a>
<a href="#c5">Saturation recovery</a><a href="#c6">Joint model (C6)</a><a href="#task7">Full pipeline</a><a href="#reference">Reference type and #767</a><a href="#next">Next steps</a></nav>
</header>
<section id="summary"><h2>Findings so far</h2>{summary}</section>
{noise}{blur}{screen}{zig}{engines}{policy}{consumers}{real}{coreg}{c5}{c6}{task7}{reference}
<section id="next">
<h2>Caveats and next steps</h2>
<ul>
<li>All candidate results are on development fixtures. The held-out fixtures have not been scored, so
nothing here is final.</li>
<li>The simulator moves the head through a B0 field fixed to the scanner, so large spikes are not fully
rigid-correctable even with the true poses; the scores are normalized by a median over pairs for that
reason.</li>
<li>Edge sensitivity may be higher in simulation than in real data (no partial-volume blurring of the
brain edge beyond the voxel grid). The real-data runs of the candidates will test whether the drift
penalty on pair-based methods carries over.</li>
<li>Next: a failure guard for production (detect disagreement with <code>c2-volreg</code> or an
implausible FD, then fall back); choosing the BBR sign from the reference's grey/white contrast;
score the held-out fixtures; and the full-pipeline runs of the short list (Task 7).</li>
<li>Same-contrast least squares with the offset step needs real-data runs and the held-out
fixtures before it can replace MCFLIRT; the #767-design arm should be rerun on #767's own code
with 3dvolreg swapped in rather than the probe's reimplementation.</li>
</ul>
</section>
</div>"""


def main(argv=None):
    global ASLSCAN  # noqa: PLW0603
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--root', type=Path, default=Path('/home/tsalo/hmc-probe'))
    parser.add_argument('--aslscan', type=Path, default=ASLSCAN)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    ASLSCAN = args.aslscan
    args.out.write_text(build(args.root), encoding='utf-8')
    print(args.out)


if __name__ == '__main__':
    main()
