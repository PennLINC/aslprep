"""Small inline-SVG charts for the report, coloured through the page's CSS tokens.

Every mark carries a ``<title>`` so hovering it shows the run and value.
"""

from html import escape

import numpy as np

#: Method -> (label, CSS token). Categorical slots in fixed order.
METHOD_STYLE = {
    'synthstrip': ('SynthStrip', 'var(--s1)'),
    'legacy': ('Legacy', 'var(--s2)'),
    'anat_pv50': ('Anatomical', 'var(--s3)'),
    'anat_nn': ('Anatomical (nearest)', 'var(--s4)'),
    'anat_pv10': ('Anatomical (≥10%)', 'var(--s5)'),
    'anat_truexfm': ('Anatomical, true transform', 'var(--s6)'),
    'hybrid': ('Dice guard', 'var(--s7)'),
    'hybrid_vol': ('Dice + volume guard', 'var(--s8)'),
    'production': ('Production run', 'var(--muted)'),
}


def _ticks(lo, hi, n=5):
    step = (hi - lo) / n
    mag = 10 ** np.floor(np.log10(step))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= step), default=step)
    start = np.ceil(lo / step) * step
    return [round(v, 10) for v in np.arange(start, hi + step / 2, step)]


def strip(rows, x_label, lo, hi, width=640, row_h=34, fmt='{:.3f}', ref=None):
    """Horizontal strip plot: one row per group, one dot per run, a median tick per row.

    ``rows``: list of (row label, colour token, [(run name, value), ...]).
    ``ref``: an x value to mark with a dashed line.
    """
    left, right, top = 190, 20, 10
    h = top + row_h * len(rows) + 44
    span = width - left - right

    def x(v):
        return left + span * (min(max(v, lo), hi) - lo) / (hi - lo)

    out = [
        f'<svg class="chart" viewBox="0 0 {width} {h}" role="img" aria-label="{escape(x_label)}">'
    ]
    for t in _ticks(lo, hi):
        out.append(
            f'<line x1="{x(t):.1f}" x2="{x(t):.1f}" y1="{top}" y2="{h - 34}" class="grid"/>'
        )
        out.append(
            f'<text x="{x(t):.1f}" y="{h - 20}" class="tick" text-anchor="middle">{fmt.format(t)}</text>'
        )
    if ref is not None:
        out.append(
            f'<line x1="{x(ref):.1f}" x2="{x(ref):.1f}" y1="{top}" y2="{h - 34}" class="refline"/>'
        )
    out.append(
        f'<text x="{left + span / 2:.1f}" y="{h - 4}" class="axis" text-anchor="middle">{escape(x_label)}</text>'
    )
    rng = np.random.default_rng(0)
    for i, (label, color, pts) in enumerate(rows):
        cy = top + row_h * i + row_h / 2
        out.append(
            f'<text x="{left - 10}" y="{cy + 4:.1f}" class="rowlab" text-anchor="end">{escape(label)}</text>'
        )
        vals = [v for _, v in pts if np.isfinite(v)]
        for name, v in pts:
            if not np.isfinite(v):
                continue
            jy = cy + rng.uniform(-row_h * 0.28, row_h * 0.28)
            clipped = ' (clipped)' if v < lo or v > hi else ''
            out.append(
                f'<circle cx="{x(v):.1f}" cy="{jy:.1f}" r="4" fill="{color}" class="dot">'
                f'<title>{escape(name)}: {fmt.format(v)}{clipped}</title></circle>'
            )
        if vals:
            m = float(np.median(vals))
            out.append(
                f'<line x1="{x(m):.1f}" x2="{x(m):.1f}" y1="{cy - row_h * 0.42:.1f}" '
                f'y2="{cy + row_h * 0.42:.1f}" class="median"><title>median {fmt.format(m)}</title></line>'
            )
    out.append('</svg>')
    return '\n'.join(out)


def bands(series, x_vals, x_label, y_label, lo, hi, width=640, height=280, fmt='{:.2f}'):
    """Median line with an interquartile band per series, over shared x values.

    ``series``: list of (label, colour token, {x: [values]}).
    """
    left, right, top, bottom = 60, 150, 12, 46
    pw, ph = width - left - right, height - top - bottom
    xs = list(x_vals)
    x0, x1 = min(xs), max(xs)

    def px(v):
        return left + pw * (v - x0) / (x1 - x0)

    def py(v):
        return top + ph * (1 - (min(max(v, lo), hi) - lo) / (hi - lo))

    out = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{escape(y_label)}">'
    ]
    for t in _ticks(lo, hi, 4):
        out.append(
            f'<line x1="{left}" x2="{left + pw}" y1="{py(t):.1f}" y2="{py(t):.1f}" class="grid"/>'
        )
        out.append(
            f'<text x="{left - 8}" y="{py(t) + 4:.1f}" class="tick" text-anchor="end">{fmt.format(t)}</text>'
        )
    out.extend(
        f'<text x="{px(v):.1f}" y="{top + ph + 18}" class="tick" text-anchor="middle">{v:g}</text>'
        for v in xs
    )
    out.append(
        f'<text x="{left + pw / 2:.1f}" y="{height - 6}" class="axis" text-anchor="middle">{escape(x_label)}</text>'
    )
    out.append(
        f'<text x="14" y="{top + ph / 2:.1f}" class="axis" text-anchor="middle" '
        f'transform="rotate(-90 14 {top + ph / 2:.1f})">{escape(y_label)}</text>'
    )
    for label, color, data in series:
        q = {v: np.percentile(data[v], [25, 50, 75]) for v in xs if len(data.get(v, []))}
        keys = sorted(q)
        poly = [f'{px(v):.1f},{py(q[v][2]):.1f}' for v in keys] + [
            f'{px(v):.1f},{py(q[v][0]):.1f}' for v in reversed(keys)
        ]
        out.append(f'<polygon points="{" ".join(poly)}" fill="{color}" opacity="0.18"/>')
        line = ' '.join(f'{px(v):.1f},{py(q[v][1]):.1f}' for v in keys)
        out.append(f'<polyline points="{line}" fill="none" stroke="{color}" stroke-width="2"/>')
        out.extend(
            f'<circle cx="{px(v):.1f}" cy="{py(q[v][1]):.1f}" r="4" fill="{color}" class="dot">'
            f'<title>{escape(label)}, {v:g}: median {fmt.format(q[v][1])} '
            f'(IQR {fmt.format(q[v][0])}–{fmt.format(q[v][2])})</title></circle>'  # noqa: RUF001
            for v in keys
        )
        last = keys[-1]
        out.append(
            f'<text x="{px(last) + 8:.1f}" y="{py(q[last][1]) + 4:.1f}" class="rowlab">{escape(label)}</text>'
        )
    out.append('</svg>')
    return '\n'.join(out)
