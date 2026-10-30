#!/usr/bin/env python3
"""Write Task 3's fixture recipes into ``scripts/coreg_mc/recipes`` from the HMC probe's.

Each new recipe is a copy of an HMC probe recipe (``scripts/hmc_probe/recipes``) with only the
settings below changed, so the copy's provenance is its source and the changes listed in its
``SOURCE.md``. The source recipes are not modified, and neither is ASLPrep's CI fixture set.

- ``raw`` variants: ``anat = 'raw'``, so ASLPrep segments the T1w itself.
- Short series with the ``rc_*`` offset: the ``mo_pcaslshort*_drift`` recipes set none.
- A second offset, with a larger rotation (8.2 degrees), on ``rc_pcasl`` and ``rc_grase``.
- Held-out still fixtures: every ``rc_*`` recipe with a new seed and a different offset, which
  the held-out fixtures share, so their pairs have a known direct registration as the
  development ones do. They are not scored until Task 5.

Generate them with ``scripts/hmc_probe/fixtures.py``'s ``generate`` and this directory as
``recipes_dir`` (``--recipes`` below does both).

Usage (in WSL)::

    python scripts/coreg_mc/make_recipes.py [--generate --data-dir <dir> [--only name ...]]
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_DIR = HERE.parent / 'hmc_probe' / 'recipes'
RECIPES_DIR = HERE / 'recipes'

#: The development fixtures' offset (tx ty tz in mm, rx ry rz in degrees).
OFFSET_DEV = [2.0, -3.0, 1.5, 3.0, -2.0, 2.5]
#: The second offset: a rotation of 8.2 degrees (norm of 6, -4, 4), for the capture range.
OFFSET_LARGE = [3.0, -4.0, 2.5, 6.0, -4.0, 4.0]
#: The held-out offset, of the development offset's size but another direction.
OFFSET_HELDOUT = [-2.5, 2.0, -1.5, -2.5, 3.0, -2.0]
SEED_HELDOUT = 20261029

RAW = ['rc_pcasl', 'rc_grase', 'st_spiral', 'st_vessels', 'st_noperf']
SHORT = ['mo_pcaslshort2_drift', 'mo_pcaslshort4_drift']
LARGE = ['rc_pcasl', 'rc_grase']
HELDOUT = ['rc_pcasl', 'rc_pcasl_flip', 'rc_pcasl_grid', 'rc_grase', 'rc_pasl_multipld']


def variants():
    """{new name: (source name, recipe settings, overlay seed or None, description)}."""
    out = {}
    for name in RAW:
        out[f'{name}_raw'] = (name, {'anat': '"raw"'}, None, "anat = 'raw' (ASLPrep segments)")

    for name in SHORT:
        out[f'{name}_off'] = (
            name,
            {'anat_offset': OFFSET_DEV},
            None,
            f'anat_offset = {OFFSET_DEV} (the rc_* offset)',
        )

    for name in LARGE:
        out[f'{name}_off2'] = (
            name,
            {'anat_offset': OFFSET_LARGE},
            None,
            f'anat_offset = {OFFSET_LARGE} (8.2 degrees)',
        )

    for name in HELDOUT:
        out[f'{name}_ho'] = (
            name,
            {'anat_offset': OFFSET_HELDOUT},
            SEED_HELDOUT,
            f'held out: anat_offset = {OFFSET_HELDOUT}, seed = {SEED_HELDOUT}',
        )
    return out


def set_key(text, key, value):
    """Set ``key = value`` in a flat TOML file, replacing the line or appending one."""
    line = f'{key} = {value}'
    new, n = re.subn(rf'(?m)^{key}\s*=.*$', line, text)
    if n > 1:
        raise ValueError(f'{key}: more than one line')
    return new if n else text.rstrip('\n') + f'\n{line}\n'


def write(name, source, settings, seed, description):
    dest = RECIPES_DIR / name
    shutil.rmtree(dest, ignore_errors=True)
    shutil.copytree(SOURCE_DIR / source, dest)
    recipe = (dest / 'recipe.toml').read_text()
    for key, value in settings.items():
        recipe = set_key(recipe, key, value)
    (dest / 'recipe.toml').write_text(recipe)
    if seed is not None:
        overlay = (dest / 'overlay.toml').read_text()
        overlay, n = re.subn(r'(?m)^seed\s*=.*$', f'seed = {seed}', overlay)
        if n != 1:
            raise ValueError(f'{source}: expected one seed line in overlay.toml')
        (dest / 'overlay.toml').write_text(overlay)

    (dest / 'SOURCE.md').write_text(
        f'# {name}\n\n'
        f'`scripts/hmc_probe/recipes/{source}`, written by `scripts/coreg_mc/make_recipes.py` '
        f'with one change: {description}.\n',
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--generate', action='store_true', help='Also generate the fixtures.')
    parser.add_argument('--data-dir', type=Path)
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args(argv)
    RECIPES_DIR.mkdir(exist_ok=True)
    table = variants()
    for name, spec in table.items():
        write(name, *spec)
        print(f'{name} <- {spec[0]}: {spec[3]}')

    if args.generate:
        sys.path.insert(0, str(HERE.parent / 'hmc_probe'))
        import fixtures

        names = args.only or sorted(table)
        fixtures.generate(names, args.data_dir, args.threads, recipes_dir=RECIPES_DIR)


if __name__ == '__main__':
    main()
