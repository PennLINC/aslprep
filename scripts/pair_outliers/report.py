# ruff: noqa: RUF001
"""Build the walkthrough HTML report from the simulation and real-data results.

Usage::

    python report.py SIM_DIR REAL_DIR OUT_HTML

Every number in the text is computed here from the JSON results.
"""

import argparse
import html
import sys
from collections import defaultdict
from pathlib import Path
from string import Template

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import simulate
import summarize as S
import svgplot as P

KINDS = list(simulate.ARTIFACTS)
KIND_LABELS = {
    'motion': 'Motion between control and label',
    'intensity': 'Slice or image intensity change',
    'labeling': 'Partial labeling failure',
    'spike': 'Excess noise',
}
SEVERITY_UNITS = {
    'motion': ('mm / degrees', [0.3, 0.6, 0.9, 1.2, 1.5]),
    'intensity': ('fractional change', None),
    'labeling': ('fraction of signal lost', [0.3, 0.475, 0.65, 0.825, 1.0]),
    'spike': ('noise SD gain', [2.0, 2.75, 3.5, 4.25, 5.0]),
}
DESIGN_ORDER = ['hcpa', 'rep4', 'rep2', 'rep1']
DESIGN_LABELS = {
    'hcpa': '2D, 5 delays × 6–15 pairs',
    'rep4': '3D, 6 delays × 4 pairs',
    'rep2': '3D, 5 delays × 2 pairs',
    'rep1': '3D, 5 delays × 1 pair',
}
#: Short protocol labels for figure panels
SHORT_DESIGN = {
    'hcpa': '2D, 6–15 per delay',
    'rep4': '3D, 4 per delay',
    'rep2': '3D, 2 per delay',
    'rep1': '3D, 1 per delay',
}
COMPARE = ['zfrac', 'zfrac_loo', 'empirical', 'model', 'hybrid']
DATASET_ORDER = ['2dmb', 'tarricone', 'lsa', 'grotzinger']
DATASET_LABELS = {
    '2dmb': '2D multiband (43 pairs, 5 delays)',
    'tarricone': '3D GRASE (12 pairs, 5 delays; 3 runs)',
    'lsa': '3D GRASE (9 pairs, 4 delays; 4 runs)',
    'grotzinger': '3D (5 pairs, 1 per delay)',
}


def pct(x, digits=0):
    return '–' if x is None or not np.isfinite(x) else f'{100 * x:.{digits}f}%'


def num(x, digits=1):
    return '–' if x is None or not np.isfinite(x) else f'{x:.{digits}f}'


def esc(text):
    return html.escape(str(text))


# ----------------------------------------------------------------------------- figures
def fig_detection(sim):
    designs = [d for d in DESIGN_ORDER if any(r['design'] == d for r in sim)]
    fig, axes = P.figure(width=7.6, height=2.8, ncols=len(designs), sharey=True)
    axes = np.atleast_1d(axes)
    width = 0.8 / len(COMPARE)
    for ax, design in zip(axes, designs, strict=True):
        runs = [r for r in sim if r['design'] == design]
        det = S.detection_by_kind(runs, COMPARE, KINDS)
        for i, method in enumerate(COMPARE):
            rates = [det[method][k][0] for k in KINDS]
            x = np.arange(len(KINDS)) + (i - (len(COMPARE) - 1) / 2) * width
            ax.bar(x, np.nan_to_num(rates), width * 0.9, color=P.METHOD_COLOR[method],
                   label=S.METHOD_LABELS[method])  # fmt: skip
        ax.set_xticks(np.arange(len(KINDS)), [k.capitalize() for k in KINDS], rotation=40,
                      ha='right', rotation_mode='anchor', fontsize=8)  # fmt: skip
        ax.set_title(SHORT_DESIGN[design], fontsize=8.5)
        ax.set_ylim(0, 1)
        ax.grid(axis='x', visible=False)
        if not any(det[m][k][0] for m in COMPARE for k in KINDS):
            ax.text(1.5, 0.5, 'No pair can be\nexcluded', ha='center', va='center',
                    fontsize=8, color=P.INK2)  # fmt: skip
    axes[0].set_ylabel('Corrupted pairs excluded')
    axes[0].yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    fig.subplots_adjust(wspace=0.12)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, -0.02),
               fontsize=8)  # fmt: skip
    return P.to_svg(fig, 'Detection rate by artifact type, design and method')


def fig_outcome(sim, outcome, ylabel, methods=('none', 'zfrac', 'zfrac_loo', 'hybrid', 'oracle')):
    designs = [d for d in DESIGN_ORDER if any(r['design'] == d for r in sim)]
    fig, axes = P.figure(width=7.6, height=2.6, ncols=len(designs))
    axes = np.atleast_1d(axes)
    for ax, design in zip(axes, designs, strict=True):
        runs = [r for r in sim if r['design'] == design]
        table = S.outcome_table(runs, list(methods), outcome)
        n_bads = sorted({cell[1] for cell in table})
        for method in methods:
            y = [table[(design, n)][method][0] for n in n_bads]
            style = {'linestyle': '--'} if method == 'oracle' else {}
            ax.plot(n_bads, y, marker='o', markersize=4, linewidth=1.6,
                    color=P.INK if method == 'oracle' else P.METHOD_COLOR[method],
                    label=S.METHOD_LABELS[method], **style)  # fmt: skip
        ax.set_xticks(n_bads)
        ax.set_title(SHORT_DESIGN[design], fontsize=8.5)
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel(ylabel)
    fig.supxlabel('Corrupted pairs in the run', fontsize=9, color=P.INK2, y=0.02)
    fig.subplots_adjust(wspace=0.32, bottom=0.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, -0.04),
               fontsize=8)  # fmt: skip
    return P.to_svg(fig, f'{ylabel} by number of corrupted pairs')


def fig_severity(sim):
    runs = [r for r in sim if r['design'] in ('hcpa', 'rep4')]
    fig, axes = P.figure(width=7.2, height=2.3, ncols=len(KINDS), sharey=True)
    for ax, kind in zip(axes, KINDS, strict=True):
        unit, edges = SEVERITY_UNITS[kind]
        if edges is None:
            values = [
                r['severity'][i]
                for r in runs
                for i in range(r['n_pairs'])
                if r['bad'][i] and r['kind'][i] == kind
            ]
            edges = list(np.quantile(values, np.linspace(0, 1, 5))) if values else [0, 1]
            edges[-1] += 1e-9
        centers = (np.array(edges[:-1]) + np.array(edges[1:])) / 2
        for method in ('zfrac', 'zfrac_loo', 'hybrid'):
            rates, _ = S.detection_by_severity(runs, method, kind, edges)
            ax.plot(centers, rates, marker='o', markersize=4, linewidth=1.6,
                    color=P.METHOD_COLOR[method], label=S.METHOD_LABELS[method])  # fmt: skip
        ax.set_title(kind.capitalize(), fontsize=8.5)
        ax.set_xlabel(unit, fontsize=8)
        ax.set_ylim(0, 1.05)
    axes[0].set_ylabel('Excluded')
    axes[0].yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    fig.subplots_adjust(wspace=0.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, bbox_to_anchor=(0.5, -0.06),
               fontsize=8)  # fmt: skip
    return P.to_svg(fig, 'Detection rate by artifact severity, hybrid versus per-delay z')


THRESHOLDS = np.arange(2.5, 8.01, 0.25)


def fig_threshold(sim):
    runs = [r for r in sim if r['design'] != 'rep1']
    fig, axes = P.figure(width=7.2, height=2.5, ncols=2)
    for method in ('model', 'hybrid'):
        clean_excl, detect = S.threshold_tradeoff(runs, method, THRESHOLDS)
        axes[0].plot(THRESHOLDS, clean_excl, color=P.METHOD_COLOR[method], linewidth=1.6,
                     label=S.METHOD_LABELS[method])  # fmt: skip
        axes[1].plot(THRESHOLDS, detect, color=P.METHOD_COLOR[method], linewidth=1.6,
                     label=S.METHOD_LABELS[method])  # fmt: skip
    for ax, title in zip(axes, ('Clean pairs above threshold', 'Corrupted pairs above threshold'),
                         strict=True):  # fmt: skip
        ax.axvline(3.5, color=P.RULE, linestyle=':', linewidth=1)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlabel('Robust z threshold')
        ax.yaxis.set_major_formatter(lambda v, _: f'{v:.0%}')
    axes[0].set_ylim(0, None)
    axes[1].set_ylim(0, 1)
    axes[1].legend(loc='lower left', fontsize=8)
    return P.to_svg(fig, 'Approximate trade-off between threshold and exclusions')


def fig_scores(example):
    """Pair scores of one run: log energies and amplitude, corrupted pairs marked."""
    scores = example['scores']['hybrid']
    n = example['n_pairs']
    bad = np.array(example['bad'])
    fig, axes = P.figure(width=7.2, height=2.6, nrows=2, sharex=True)
    x = np.arange(n)
    for ax, key, label in (
        (axes[0], 'energy_smooth', 'log smoothed energy'),
        (axes[1], 'amplitude', 'amplitude'),
    ):
        values = np.log(scores[key]) if key.startswith('energy') else np.array(scores[key])
        excluded = ~np.array(example['keep']['hybrid'])
        ax.vlines(x, np.median(values), values, color=P.GRID, linewidth=1)
        ax.scatter(x[~bad], values[~bad], s=14, color=P.METHOD_COLOR['none'], zorder=3,
                   label='Clean pair')  # fmt: skip
        ax.scatter(x[bad], values[bad], s=28, color=P.METHOD_COLOR['zfrac'], zorder=4,
                   marker='D', label='Corrupted pair')  # fmt: skip
        ax.scatter(x[excluded], values[excluded], s=80, facecolors='none', zorder=5,
                   edgecolors=P.METHOD_COLOR['hybrid'], linewidths=1.4,
                   label='Excluded by hybrid')  # fmt: skip
        ax.set_ylabel(label, fontsize=8)
    plds = example['pld']
    ticks = [i for i in range(n) if i == 0 or plds[i] != plds[i - 1]]
    axes[1].set_xticks(ticks, [f'{plds[i]:g} s' for i in ticks])
    axes[1].set_xlabel('Pairs, in acquisition order (labels mark the first pair at each delay)')
    axes[0].legend(loc='lower left', bbox_to_anchor=(0, 1.02), ncol=3, fontsize=7.5,
                   borderaxespad=0)  # fmt: skip
    return P.to_svg(fig, 'Scores of every pair in one simulated run')


# ----------------------------------------------------------------------------- tables
def outcome_rows(sim, methods):
    rows = []
    for design in DESIGN_ORDER:
        runs = [r for r in sim if r['design'] == design]
        if not runs:
            continue
        tables = {o: S.outcome_table(runs, methods, o) for o in S.OUTCOMES}
        rows.extend(
            (design, cell[1], {o: tables[o][cell] for o in S.OUTCOMES})
            for cell in sorted(tables['dm_rmse_pct'])
        )
    return rows


def table_html(headers, rows, cls=''):
    head = ''.join(f'<th>{h}</th>' for h in headers)
    body = ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in row) + '</tr>' for row in rows)
    return (
        f'<div class="table-wrap"><table class="{cls}"><thead><tr>{head}</tr></thead>'
        f'<tbody>{body}</tbody></table></div>'
    )


# ----------------------------------------------------------------------------- page
def _real_detection(real, method, dataset=None):
    """Share of injected pairs excluded, and of newly excluded clean pairs, on real data."""
    observed = {r['run']: r['keep'][method] for r in real if r['n_bad'] == 0}
    hits, fp, n_clean = [], 0, 0
    for r in real:
        if r['n_bad'] == 0 or r['run'] not in observed:
            continue
        if dataset is not None and r['run'].split('_')[0] != dataset:
            continue
        keep, bad = np.array(r['keep'][method]), np.array(r['bad'])
        base = np.array(observed[r['run']])
        hits += list(~keep[bad])
        fp += int((~keep & ~bad & base).sum())
        n_clean += int((~bad & base).sum())
    return (float(np.mean(hits)) if hits else np.nan), fp / max(n_clean, 1)


def conclusions(sim, real, rates):
    """The closing section, with its numbers computed from the results."""
    design_rates = {
        d: S.exclusion_rates([r for r in sim if r['design'] == d], ['hybrid', 'zfrac_loo'])
        for d in DESIGN_ORDER
    }
    clean = {
        d: {
            m: np.mean(
                [
                    r['metrics'][m]['dm_rmse_pct']
                    for r in sim
                    if r['design'] == d and r['n_bad'] == 0
                ]
            )
            for m in ('none', 'hybrid')
        }
        for d in ('hcpa', 'rep4')
    }
    clean_cost = {d: 100 * (v['hybrid'] / v['none'] - 1) for d, v in clean.items()}
    real_text = ''
    if real:

        def mean_metric(dataset, n_bad, method, key):
            runs = [r for r in real if r['run'].startswith(dataset) and r['n_bad'] == n_bad]
            return np.nanmean([r['metrics'][method][key] for r in runs])

        mb_det, _ = _real_detection(real, 'hybrid', '2dmb')
        mb_zl, _ = _real_detection(real, 'zfrac_loo', '2dmb')
        mb_heavy = max(r['n_bad'] for r in real if r['run'].startswith('2dmb'))
        three_d = [
            _real_detection(real, 'hybrid', d) for d in ('tarricone', 'lsa')
        ]  # detection, false exclusion
        real_text = f"""
<li>On the real 2D multiband run with injected artifacts, the hybrid detector excluded
{pct(mb_det)} of the injected pairs (leave-one-out z rule: {pct(mb_zl)}). With {mb_heavy}
corrupted pairs, it cut the mean ΔM error from
{num(mean_metric('2dmb', mb_heavy, 'none', 'dm_rmse_pct'), 0)}% to
{num(mean_metric('2dmb', mb_heavy, 'hybrid', 'dm_rmse_pct'), 0)}% and the median regional CBF
error from {num(mean_metric('2dmb', mb_heavy, 'none', 'region_cbf_mae_pct'), 1)}% to
{num(mean_metric('2dmb', mb_heavy, 'hybrid', 'region_cbf_mae_pct'), 1)}% (oracle:
{num(mean_metric('2dmb', mb_heavy, 'oracle', 'dm_rmse_pct'), 0)}% and
{num(mean_metric('2dmb', mb_heavy, 'oracle', 'region_cbf_mae_pct'), 1)}%).</li>
<li>On the real 3D GRASE runs with two or three pairs per delay, it excluded only
{pct(three_d[0][0])} and {pct(three_d[1][0])} of the injected pairs, newly excluded
{pct(three_d[0][1], 1)} and {pct(three_d[1][1], 1)} of clean pairs, and did not improve the
fits. The simulations with two pairs per delay were more optimistic, so the simulated 3D noise
or artifact severities do not match these data; at this depth of sampling the method should be
treated as a QC score rather than an exclusion rule.</li>
"""
    return f"""
<h2>Conclusions</h2>
<ul>
<li>Model-based pair exclusion works for multi-delay ASL when each delay has several repeats.
With four or more pairs per delay, the hybrid detector excluded
{pct(design_rates['rep4']['hybrid']['detect'])} (3D, four per delay) and
{pct(design_rates['hcpa']['hybrid']['detect'])} (2D, HCP-Aging-like) of corrupted pairs, and its
ΔM and kinetic fits came close to those of the oracle that knows which pairs are corrupted.</li>
<li>It finds more corrupted pairs than either model-free per-delay rule
({pct(rates['hybrid']['detect'])} overall, against {pct(rates['zfrac']['detect'])} and
{pct(rates['zfrac_loo']['detect'])}). The kinetic model helps most where a delay has few pairs to
compare: the published per-delay rule cannot act at all with four or fewer repeats.</li>
<li>With two pairs per delay it still finds {pct(design_rates['rep2']['hybrid']['detect'])} of
corrupted pairs, but only by leaning on the model, because two pairs cannot outvote each other.
With one pair per delay nothing can be excluded without losing a delay, so the method only
scores the pairs.</li>
<li>The cost is a small rate of false exclusions: {pct(rates['hybrid']['clean_excl'], 1)} of clean pairs
overall. In runs without corruption this raised the mean ΔM error by
{num(clean_cost['hcpa'], 1)}% (2D) and {num(clean_cost['rep4'], 1)}% (3D, four per delay)
relative to keeping every pair.</li>
{real_text}
</ul>

<h3>Recommendation</h3>
<p>Add this to ASLPrep for multi-delay data, between delta-M extraction and the kinetic fit in
<code>ComputeCBF</code>, and always write the per-pair scores and decisions to a TSV (the
observation-level table proposed in PennLINC/aslprep#683), so that exclusions can be audited
and the scores are available as QC even when nothing is excluded. Make exclusion opt-in at
first, and recommend it for protocols with four or more pairs per delay, where both the
simulations and the real 2D run support it; with fewer, report the scores and keep every pair
unless the user asks otherwise. Start from the defaults tested here: κ = 2, threshold 3.5,
smoothing σ = 4 mm, and at least one pair kept per delay. Check how many other pairs the
same-delay correction should need (three here; two may be better, see above) on independent
data first, such as the aslscan motion fixtures or the real runs of the
<code>asl-outlier-eval</code> study.</p>

<h3>Limitations</h3>
<ul>
<li>The kinetic-model part of each pair's prediction is not leave-one-out. A corrupted pair
pulls the fit toward itself, most at delays with few pairs; the robust first fit and the
iteration limit this, but do not remove it.</li>
<li>The thresholds, the smoothing width and the three-pair rule were set on these simulations,
and the three-pair rule was chosen after the first run (see above).</li>
<li>The simulations are PCASL at 4 mm, with a forward model from the same family as the fitted
one (plus dispersion, tissue T1 and physiological noise). PASL is covered only by the unit tests,
and time-encoded (Hadamard) acquisitions not at all.</li>
<li>The real-data references are the runs as acquired, so a pair that was already corrupted
counts as clean, and the severity of injected noise spikes is nominal (scaled by a rough noise
estimate).</li>
<li>Exclusion is per pair. 2D slice dropouts are handled by excluding the whole pair; excluding
single slices would need per-voxel weights in the fit.</li>
<li>Robust or heavy-tailed likelihoods inside a Bayesian fit, an alternative closer to BASIL,
were not tested.</li>
</ul>
"""


def build(sim, real, out, sim_v1=(), prep_dir=None):
    all_methods = ['none', 'zfrac', 'zfrac_loo', 'empirical', 'model', 'hybrid', 'hybrid_raw',
                   'hybrid_noamp', 'oracle']  # fmt: skip
    rates = S.exclusion_rates(sim, all_methods)
    by_design = {
        d: S.exclusion_rates([r for r in sim if r['design'] == d], all_methods)
        for d in DESIGN_ORDER
    }
    clean_runs = S.clean_run_exclusions(sim, all_methods)
    n_runs = len(sim)
    n_seeds = len({r['seed'] for r in sim})

    # Headline numbers
    hyb, zf, zl = rates['hybrid'], rates['zfrac'], rates['zfrac_loo']
    heavy = [r for r in sim if r['design'] == 'hcpa' and r['n_bad'] == max(
        rr['n_bad'] for rr in sim if rr['design'] == 'hcpa')]  # fmt: skip
    dm_heavy = {m: np.mean([r['metrics'][m]['dm_rmse_pct'] for r in heavy]) for m in
                ('none', 'zfrac', 'zfrac_loo', 'hybrid', 'oracle')}  # fmt: skip
    clean_hcpa = [r for r in sim if r['design'] == 'hcpa' and r['n_bad'] == 0]
    dm_clean = {m: np.mean([r['metrics'][m]['dm_rmse_pct'] for r in clean_hcpa]) for m in
                ('none', 'hybrid')}  # fmt: skip

    # Detection table, all designs pooled
    det = S.detection_by_kind(sim, all_methods[1:-1], KINDS)
    det_rows = [
        [S.METHOD_LABELS[method]]
        + [pct(det[method][k][0]) for k in KINDS]
        + [
            pct(rates[method]['detect']),
            pct(rates[method]['clean_excl'], 1),
            pct(clean_runs[method]),
        ]
        for method in all_methods[1:-1]
    ]
    det_table = table_html(
        ['Method'] + [k.capitalize() for k in KINDS]
        + ['All artifacts', 'Clean pairs excluded', 'Clean runs losing a pair'],
        det_rows, 'num',
    )  # fmt: skip
    n_kind = {k: det['hybrid'][k][1] for k in KINDS}

    # False exclusions by design
    clean_excl_methods = ('zfrac', 'zfrac_loo', 'empirical', 'model', 'hybrid')
    clean_excl_rows = [
        [DESIGN_LABELS[design]]
        + [pct(by_design[design][m]['clean_excl'], 1) for m in clean_excl_methods]
        + [pct(by_design[design][m]['detect']) for m in ('zfrac_loo', 'hybrid')]
        for design in DESIGN_ORDER
    ]
    clean_excl_table = table_html(
        ['Design', 'Per-delay z', 'Per-delay z (LOO)', 'Same-delay median', 'Model only',
         'Hybrid', 'Detection: z (LOO)', 'Detection: hybrid'], clean_excl_rows, 'num')  # fmt: skip

    # Outcome tables
    out_methods = ['none', 'zfrac', 'zfrac_loo', 'hybrid', 'oracle']
    orows = []
    for design, n_bad, cells in outcome_rows(sim, out_methods):
        row = [DESIGN_LABELS[design], str(n_bad)]
        for outcome, digits in (('dm_rmse_pct', 0), ('region_cbf_rmse', 1),
                                ('region_att_rmse', 2)):  # fmt: skip
            vals = cells[outcome]
            best = min(vals[m][0] for m in out_methods if m != 'oracle')
            for m in out_methods:
                v = vals[m][0]
                text = num(v, digits)
                if m != 'oracle' and v == best:
                    text = f'<strong>{text}</strong>'
                row.append(text)
        orows.append(row)
    sub = ['None', 'Z', 'Z LOO', 'Hybrid', 'Oracle']
    outcome_table = table_html(
        ['Design', 'Bad pairs']
        + [f'ΔM err. %<br>{s}' for s in sub]
        + [f'Reg. CBF<br>{s}' for s in sub]
        + [f'Reg. ATT (s)<br>{s}' for s in sub],
        orows,
        'num dense',
    )

    # Ablations: mean outcome differences against hybrid, pooled over contaminated runs
    contaminated = [r for r in sim if r['n_bad'] > 0]
    abl_rows = []
    for method in ('empirical', 'model', 'hybrid_raw', 'hybrid_noamp', 'hybrid'):
        dm = np.mean([r['metrics'][method]['dm_rmse_pct'] for r in contaminated])
        cbf = np.mean([r['metrics'][method]['region_cbf_rmse'] for r in contaminated])
        abl_rows.append([S.METHOD_LABELS[method], pct(rates[method]['detect']),
                         pct(rates[method]['clean_excl'], 1), num(dm, 1), num(cbf, 1)])  # fmt: skip
    abl_table = table_html(
        ['Variant', 'Detection', 'Clean pairs excluded', 'ΔM error %', 'Regional CBF RMSE'],
        abl_rows,
        'num',
    )

    # What the amplitude test costs and buys
    amp_note = (
        f'The amplitude test is what catches partial labeling failures '
        f'({pct(det["hybrid"]["labeling"][0])} with it, '
        f'{pct(det["hybrid_noamp"]["labeling"][0])} without), but it also causes most of the '
        f'false exclusions in the 3D protocol with four pairs per delay '
        f'({pct(by_design["rep4"]["hybrid"]["clean_excl"], 1)} of clean pairs with it, '
        f'{pct(by_design["rep4"]["hybrid_noamp"]["clean_excl"], 1)} without). Smoothed energy is '
        f'what catches intensity changes ({pct(det["hybrid"]["intensity"][0])} with it, '
        f'{pct(det["hybrid_raw"]["intensity"][0])} without).'
    )

    # Example run for the score figure: heavy HCP-A run with the most artifact kinds
    example = max(heavy, key=lambda r: len({k for k in r['kind'] if k}))

    if real:
        real_html = real_section(real, prep_dir)
    else:
        real_html = '<p>Real-data results are not available.</p>'

    # First run (same-delay correction from one other pair) against this run
    version_stats = {}
    version_rows = []
    for design in ('rep4', 'rep2'):
        for label, results in (('First run', sim_v1), ('This run', sim)):
            runs = [r for r in results if r['design'] == design]
            if not runs:
                continue
            rate = S.exclusion_rates(runs, ['hybrid'])['hybrid']
            bad_runs = [r for r in runs if r['n_bad'] > 0]
            version_stats[(design, label)] = (
                rate['detect'],
                rate['clean_excl'],
                np.mean([r['metrics']['hybrid']['region_cbf_rmse'] for r in bad_runs]),
            )
            version_rows.append([
                DESIGN_LABELS[design], label, pct(rate['detect']), pct(rate['clean_excl'], 1),
                num(np.mean([r['metrics']['hybrid']['dm_rmse_pct'] for r in bad_runs]), 1),
                num(np.mean([r['metrics']['none']['dm_rmse_pct'] for r in bad_runs]), 1),
                num(np.mean([r['metrics']['hybrid']['region_cbf_rmse'] for r in bad_runs]), 1),
                num(np.mean([r['metrics']['none']['region_cbf_rmse'] for r in bad_runs]), 1),
            ])  # fmt: skip
    version_note = ''
    if ('rep4', 'First run') in version_stats:

        def change(design, i, fmt):
            return f'{fmt(version_stats[(design, "First run")][i])} to ' + fmt(
                version_stats[(design, 'This run')][i]
            )

        version_note = (
            'The rule did what it was meant to with two pairs per delay (detection from '
            f'{change("rep2", 0, pct)}, clean pairs excluded from '
            f'{change("rep2", 1, lambda x: pct(x, 1))}, regional CBF RMSE from '
            f'{change("rep2", 2, num)}), but it hurt with four (detection from '
            f'{change("rep4", 0, pct)}, clean pairs excluded from '
            f'{change("rep4", 1, lambda x: pct(x, 1))}, regional CBF RMSE from '
            f'{change("rep4", 2, num)}). With four pairs per delay, once one pair at a delay is '
            'excluded only two others remain, so the correction switches off and those pairs '
            'are predicted from the model alone. Requiring two other pairs instead of three '
            "would probably keep the first run's behaviour with four pairs per delay and this "
            "run's with two, but choosing it from these same simulations would be a third "
            'change tuned on them. It should be checked on independent data before it is '
            'adopted.'
        )
    version_table = table_html(
        ['Design', 'Run', 'Detection', 'Clean pairs excluded', 'ΔM error % (hybrid)',
         'ΔM error % (none)', 'Reg. CBF RMSE (hybrid)', 'Reg. CBF RMSE (none)'],
        version_rows, 'num')  # fmt: skip

    page = TEMPLATE.substitute(
        n_runs=n_runs,
        n_seeds=n_seeds,
        hyb_detect=pct(hyb['detect']),
        zf_detect=pct(zf['detect']),
        zl_detect=pct(zl['detect']),
        hyb_clean_excl=pct(hyb['clean_excl'], 1),
        zf_clean_excl=pct(zf['clean_excl'], 1),
        zl_clean_excl=pct(zl['clean_excl'], 1),
        dm_none=num(dm_heavy['none'], 0),
        dm_zfrac=num(dm_heavy['zfrac'], 0),
        dm_zloo=num(dm_heavy['zfrac_loo'], 0),
        dm_hybrid=num(dm_heavy['hybrid'], 0),
        dm_oracle=num(dm_heavy['oracle'], 0),
        dm_clean_none=num(dm_clean['none'], 1),
        dm_clean_hybrid=num(dm_clean['hybrid'], 1),
        n_heavy=heavy[0]['n_bad'],
        n_kind=', '.join(f'{n_kind[k]} {k}' for k in KINDS),
        det_table=det_table,
        clean_excl_table=clean_excl_table,
        outcome_table=outcome_table,
        abl_table=abl_table,
        fig_detection=fig_detection(sim),
        fig_dm=fig_outcome(sim, 'dm_rmse_pct', 'Mean ΔM error (%)'),
        fig_cbf=fig_outcome(sim, 'region_cbf_rmse', 'Regional CBF RMSE'),
        fig_severity=fig_severity(sim),
        fig_threshold=fig_threshold(sim),
        fig_scores=fig_scores(example),
        example_desc=esc(
            ', '.join(
                f'pair {i + 1} ({example["kind"][i]}, {example["pld"][i]:g} s)'
                for i in range(example['n_pairs'])
                if example['bad'][i]
            )
        ),
        real_section=real_html,
        version_table=version_table,
        version_note=version_note,
        amp_note=amp_note,
        empirical_detect=pct(rates['empirical']['detect']),
        conclusions=conclusions(sim, real, rates),
    )
    Path(out).write_text(page, encoding='utf-8')


def real_section(real, prep_dir=None):
    observed = sorted([r for r in real if r['n_bad'] == 0], key=lambda r: r['run'])
    injected = [r for r in real if r['n_bad'] > 0]
    methods = ['zfrac', 'zfrac_loo', 'model', 'hybrid']

    obs_rows = []
    for r in observed:
        plds = np.array(r['pld'])
        uniq, counts = np.unique(plds, return_counts=True)
        design = ', '.join(f'{c}×{p:g}' for p, c in zip(uniq, counts, strict=True))

        def excluded(method, r=r, plds=plds):
            idx = np.flatnonzero(~np.array(r['keep'][method]))
            return ', '.join(f'{i + 1} ({plds[i]:g} s)' for i in idx) or 'none'

        obs_rows.append([esc(r['run']), esc(design), excluded('zfrac'), excluded('zfrac_loo'),
                         excluded('hybrid'), num(r['metrics']['hybrid']['dm_rmse_pct'], 1),
                         num(r['metrics']['hybrid']['region_cbf_mae_pct'], 1)])  # fmt: skip
    obs_table = table_html(
        ['Run', 'Pairs × delay (s)', 'Excluded: per-delay z', 'Excluded: z (LOO)',
         'Excluded: hybrid', 'ΔM change, hybrid (%)', 'Regional CBF change, hybrid (%)'],
        obs_rows)  # fmt: skip

    # Injection: detection by kind and false exclusions, relative to the run as acquired
    obs_keep = {r['run']: r['keep'] for r in observed}
    datasets = sorted({r['run'].split('_')[0] for r in injected}, key=DATASET_ORDER.index)
    inj_rows = []
    for dataset in datasets:
        for method in methods:
            hits = defaultdict(list)
            fp = n_clean = 0
            for r in injected:
                if r['run'].split('_')[0] != dataset:
                    continue
                keep = np.array(r['keep'][method])
                bad = np.array(r['bad'])
                base = np.array(obs_keep[r['run']][method])
                for i in np.flatnonzero(bad):
                    hits[r['kind'][i]].append(not keep[i])
                # Clean pairs excluded now but kept when nothing was injected
                fp += int((~keep & ~bad & base).sum())
                n_clean += int((~bad & base).sum())
            all_hits = [h for k in hits for h in hits[k]]
            inj_rows.append(
                [esc(DATASET_LABELS[dataset]), S.METHOD_LABELS[method]]
                + [pct(np.mean(hits[k])) if hits[k] else '–' for k in KINDS]
                + [pct(np.mean(all_hits)), pct(fp / max(n_clean, 1), 1)]
            )
    inj_table = table_html(
        ['Data', 'Method'] + [k.capitalize() for k in KINDS]
        + ['All', 'Newly excluded clean pairs'], inj_rows, 'num')  # fmt: skip

    # Injection outcomes by dataset
    groups = defaultdict(list)
    for r in injected:
        groups[(DATASET_ORDER.index(r['run'].split('_')[0]), r['n_bad'])].append(r)
    out_rows = []
    for (i_dataset, n_bad), runs in sorted(groups.items()):
        row = [esc(DATASET_LABELS[DATASET_ORDER[i_dataset]]), str(n_bad), str(len(runs))]
        row += [
            num(np.nanmean([r['metrics'][m][outcome] for r in runs]), digits)
            for outcome, digits in (('dm_rmse_pct', 0), ('region_cbf_mae_pct', 1))
            for m in ('none', 'zfrac', 'zfrac_loo', 'hybrid', 'oracle')
        ]
        out_rows.append(row)
    sub = ['None', 'Z', 'Z LOO', 'Hybrid', 'Oracle']
    out_table = table_html(
        ['Dataset', 'Injected', 'Runs'] + [f'ΔM err. %<br>{s}' for s in sub]
        + [f'Reg. CBF median err. %<br>{s}' for s in sub], out_rows, 'num dense')  # fmt: skip

    # Do the hybrid detector's exclusions on the acquired runs coincide with motion?
    notes = []
    for r in observed:
        excluded = np.flatnonzero(~np.array(r['keep']['hybrid']))
        motion_file = prep_dir / f'{r["run"]}_motion.par' if prep_dir else None
        if not excluded.size or motion_file is None or not motion_file.exists():
            continue
        par = np.loadtxt(motion_file)
        # Displacement between each pair's control and label volumes (rotations at 50 mm)
        moved = np.hstack([par[:, :3] * 50, par[:, 3:]])
        within = np.abs(moved[1::2] - moved[0::2]).sum(axis=1)
        severity = np.array(r['scores']['hybrid']['severity'])
        for i in excluded:
            margin = 'just above' if severity[i] < 4.5 else 'above'
            cause = (
                'so the exclusion does not reflect motion'
                if within[i] <= np.median(within)
                else 'more than most pairs in the run'
            )
            notes.append(
                f'In {esc(r["run"])}, the hybrid detector excluded pair {i + 1} '
                f'({r["pld"][i]:g} s) with a severity of {severity[i]:.1f}, {margin} the '
                f'threshold of 3.5. Its control and label volumes moved {within[i]:.2f} mm '
                f'relative to each other, against a median of {np.median(within):.2f} mm for '
                f'the run, {cause}.'
            )
    observed_note = f'<p>{" ".join(notes)}</p>' if notes else ''

    n_inj = len(injected)
    return REAL_TEMPLATE.substitute(
        n_runs=len(observed), n_inj=n_inj, obs_table=obs_table, inj_table=inj_table,
        out_table=out_table, observed_note=observed_note,
    )  # fmt: skip


TEMPLATE = Template((Path(__file__).parent / 'report_template.html').read_text(encoding='utf-8'))
REAL_TEMPLATE = Template((Path(__file__).parent / 'report_real.html').read_text(encoding='utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('sim_dir', type=Path)
    parser.add_argument('real_dir', type=Path)
    parser.add_argument('out', type=Path)
    parser.add_argument('--sim-v1', type=Path, help='results of the first simulation run')
    parser.add_argument('--prep-dir', type=Path, help='real runs prepared by real_prep.py')
    args = parser.parse_args()
    sim = S.load(args.sim_dir)
    real = S.load(args.real_dir) if args.real_dir.exists() else []
    sim_v1 = S.load(args.sim_v1) if args.sim_v1 else []
    build(sim, real, args.out, sim_v1, args.prep_dir)
    print(f'wrote {args.out} from {len(sim)} simulated and {len(real)} real results')


if __name__ == '__main__':
    main()
