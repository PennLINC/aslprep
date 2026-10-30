"""Matplotlib figures as inline SVG that follow the report page's light and dark tokens.

Figures are drawn with placeholder colors, which are replaced by CSS custom properties after
rendering, so axes, text and series take the page's theme.
"""

import io
import re

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

#: Placeholder color -> CSS token. The placeholders are unlikely to be drawn by accident.
TOKENS = {
    '#010101': 'var(--ink)',
    '#020202': 'var(--ink-2)',
    '#030303': 'var(--rule)',
    '#040404': 'var(--grid)',
    '#050505': 'var(--surface)',
    '#0a0001': 'var(--m-none)',
    '#0a0002': 'var(--m-zfrac)',
    '#0a0003': 'var(--m-hybrid)',
    '#0a0004': 'var(--m-oracle)',
    '#0a0005': 'var(--m-model)',
    '#0a0006': 'var(--m-empirical)',
    '#0a0007': 'var(--m-zfrac-loo)',
}

INK, INK2, RULE, GRID, SURFACE = '#010101', '#020202', '#030303', '#040404', '#050505'

#: Series placeholder of each method
METHOD_COLOR = {
    'none': '#0a0001',
    'zfrac': '#0a0002',
    'hybrid': '#0a0003',
    'oracle': '#0a0004',
    'model': '#0a0005',
    'empirical': '#0a0006',
    'zfrac_loo': '#0a0007',
}


def figure(width=6.4, height=3.2, **kwargs):
    """Make a figure styled for the report."""
    plt.rcParams.update(
        {
            'font.family': 'sans-serif',
            'font.size': 9,
            'axes.edgecolor': RULE,
            'axes.labelcolor': INK2,
            'axes.titlecolor': INK,
            'axes.titlesize': 9.5,
            'axes.titleweight': 'bold',
            'axes.titlelocation': 'left',
            'xtick.color': INK2,
            'ytick.color': INK2,
            'text.color': INK,
            'axes.grid': True,
            'grid.color': GRID,
            'grid.linewidth': 0.6,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'legend.frameon': False,
            'svg.fonttype': 'none',
        }
    )
    fig, axes = plt.subplots(figsize=(width, height), **kwargs)
    fig.patch.set_alpha(0)
    for ax in fig.axes:
        ax.set_facecolor('none')
        ax.set_axisbelow(True)
    return fig, axes


def to_svg(fig, label):
    """Render a figure to inline SVG, with tokens in place of the placeholder colors."""
    buffer = io.StringIO()
    fig.savefig(buffer, format='svg', bbox_inches='tight', transparent=True)
    plt.close(fig)
    svg = buffer.getvalue()
    svg = svg[svg.index('<svg') :]
    for placeholder, token in TOKENS.items():
        svg = re.sub(placeholder, token, svg, flags=re.IGNORECASE)
    # Scale to the container, and label the figure for screen readers.
    svg = re.sub(r'width="[\d.]+pt" height="[\d.]+pt"', '', svg, count=1)
    return svg.replace('<svg ', f'<svg role="img" aria-label="{label}" class="fig" ', 1)
