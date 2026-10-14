#!/usr/bin/env python3
"""Extend the real-data panel with 20 more subjects per family (plan: "Extending the panel").

For each extended dataset, eligible subjects are listed from the DataLad clone that the zoo's
own scripts made (``~/.cache/aslprep-zoo/<dataset>``), 20 are drawn with a fixed seed (the zoo's
exemplar subject excluded), only the files ASLPrep needs are fetched, and a BIDS dataset is
built under ``--out`` with the same fixes the zoo's curation applies to that dataset (and the
probe's adaptations from ``panel_config.py``). The fetched content is dropped from the clone
afterwards, to save space. The zoo itself is never touched.

Each dataset gets an ``EXTENDED.md`` listing the subjects and every change, and the
development / held-out split is written to ``<out>/split.json``.

Usage::

    python scripts/hmc_probe/extend_panel.py --out <dir> [--only ds007116 ...] [--n 20]
"""

import argparse
import json
import random
import re
import shutil
import subprocess
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

CACHE = Path.home() / '.cache' / 'aslprep-zoo'
SEED = 20261014
N_HELD_OUT = 5


def _cmd_datalad(ds_dir, op, paths):
    """Run ``datalad get`` or ``drop`` in the ``dlad`` environment, in batches."""
    for i in range(0, len(paths), 200):
        batch = [str(p) for p in paths[i : i + 200]]
        cmd = ['micromamba', 'run', '-n', 'dlad', 'datalad', op, '-d', str(ds_dir)]
        # The OpenNeuro S3 copy cannot be verified from here (numcopies is 2), but it is where
        # the content was fetched from and can be fetched again.
        cmd += (['-J', '4'] if op == 'get' else ['--reckless', 'availability']) + batch
        subprocess.run(cmd, check=True, capture_output=True)


def _read_json(path):
    return json.loads(Path(path).read_text())


def _write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2) + '\n')


# ---------------------------------------------------------------------------------------------
# Per-dataset rules: which runs, which anatomical image, which metadata fixes
# ---------------------------------------------------------------------------------------------
class Spec:
    """How one dataset is extended. Subclasses override what differs."""

    name = ''
    exemplar = ''
    #: Regular expression on the ASL file name (without subject/session) of eligible runs.
    asl_pattern = r'(.*_)?asl\.nii\.gz'
    #: Regular expression on anatomical file names: the T1w to keep (first match by name).
    t1w_pattern = r'_T1w\.nii\.gz$'
    #: Minimum number of sessions with an eligible run.
    min_sessions = 1
    #: Use at most this many sessions per subject (the first ones by name).
    max_sessions = None
    separate_m0 = True
    root_aslcontext = False
    notes = ()

    def __init__(self):
        self.src = CACHE / self.name

    def runs(self, sub_dir):
        """Eligible ASL runs of a subject, grouped by session directory."""
        out = {}
        for asl in sorted(sub_dir.glob('**/perf/*_asl.nii.gz')):
            rest = re.sub(r'^sub-[^_]+_(ses-[^_]+_)?', '', asl.name)
            if re.fullmatch(self.asl_pattern, rest):
                out.setdefault(asl.parent.parent, []).append(asl)
        sessions = sorted(out)
        if self.max_sessions:
            sessions = sessions[: self.max_sessions]
        return {s: out[s] for s in sessions}

    def t1w(self, ses_dir):
        anat = ses_dir / 'anat'
        cands = sorted(p for p in anat.glob('*.nii.gz') if re.search(self.t1w_pattern, p.name))
        return cands[0] if cands else None

    def t1w_name(self, src):
        return src.name

    @staticmethod
    def _stems(asl):
        stem = asl.name.replace('_asl.nii.gz', '')
        return (stem, re.sub(r'_echo-\d+', '', stem), re.sub(r'_(acq|run|dir)-[^_]+', '', stem))

    def m0(self, asl):
        if not self.separate_m0:
            return None
        for stem in self._stems(asl):
            cand = asl.with_name(f'{stem}_m0scan.nii.gz')
            if cand.exists() or cand.is_symlink():
                return cand
        return None

    def aslcontext(self, asl):
        if self.root_aslcontext:
            return self.src / 'aslcontext.tsv'
        for stem in self._stems(asl):
            p = asl.with_name(f'{stem}_aslcontext.tsv')
            if p.exists():
                return p
        return None

    def eligible(self, sub_dir):
        runs = self.runs(sub_dir)
        if len(runs) < self.min_sessions:
            return False
        for ses_dir, asls in runs.items():
            if self.t1w(ses_dir) is None:
                return False
            for asl in asls:
                if self.aslcontext(asl) is None:
                    return False
                if self.separate_m0 and self.m0(asl) is None:
                    return False
        return True

    def fix_asl_meta(self, meta, asl):
        return []

    def fix_context(self, context, asl, data):
        return context, data, []


class DS008719(Spec):
    name, exemplar = 'ds008719', 'sub-050'
    asl_pattern = r'echo-1_asl\.nii\.gz'
    notes = ('Echo 1 only: ASLPrep outputs drop the echo entity, so two echoes would collide.',)


class DS004856(Spec):
    name, exemplar = 'ds004856', 'sub-3131'
    t1w_pattern = r'_acq-MPRAGE_run-1_T1w\.nii\.gz$'
    separate_m0 = False
    max_sessions = 1

    def runs(self, sub_dir):
        runs = super().runs(sub_dir)
        return {s: a for s, a in runs.items() if s.name == 'ses-wave1'}

    def aslcontext(self, asl):
        p = asl.parent / f'{asl.parent.parent.parent.name}_{asl.parent.parent.name}_aslcontext.tsv'
        return p if p.exists() else None

    def fix_asl_meta(self, meta, asl):
        tr = meta['RepetitionTime']
        pair = tr > 6
        meta['RepetitionTimePreparation'] = round(tr / (2 if pair else 1), 5)
        return [f'RepetitionTimePreparation -> {meta["RepetitionTimePreparation"]} s']

    def fix_context(self, context, asl, data):
        """Choose the shipped alternating order or the stacked order by perfusion contrast."""
        n = len(context)
        mean = data.mean(-1)
        brain = mean > np.percentile(mean, 70)

        def contrast(ctx):
            c = [i for i, t in enumerate(ctx) if t == 'control']
            lab = [i for i, t in enumerate(ctx) if t == 'label']
            d = data[..., c].mean(-1) - data[..., lab].mean(-1)
            v = d[brain]
            return v.mean() / v.std()

        stacked = ['control'] * (n // 2) + ['label'] * (n // 2)
        shipped, alt = contrast(context), contrast(stacked)
        if alt > max(shipped, 0.2) and shipped < 0.2:
            return (
                stacked,
                data,
                [f'aslcontext: stacked order (contrast {alt:.2f} vs {shipped:.2f})'],
            )
        return (
            context,
            data,
            [f'aslcontext as shipped (contrast {shipped:.2f}; stacked {alt:.2f})'],
        )


class DS004146(Spec):
    name, exemplar = 'ds004146', 'sub-0337'
    t1w_pattern = r'_UNIT1_denoised\.nii\.gz$'
    min_sessions = 2
    separate_m0 = False
    notes = (
        'T1w: a copy of the denoised MP2RAGE UNIT1, as in the zoo.',
        'The separate M0 is not copied: M0Type is Included and ASLPrep raises an error when a '
        'separate M0 is also associated.',
    )

    def t1w_name(self, src):
        return src.name.replace('_UNIT1_denoised.nii.gz', '_rec-denoised_T1w.nii.gz')


class DS006148(Spec):
    name, exemplar = 'ds006148', 'sub-C00045'
    asl_pattern = r'dir-AP_asl\.nii\.gz'
    t1w_pattern = r'_ses-S001_T1w\.nii\.gz$'
    separate_m0 = False
    root_aslcontext = True
    max_sessions = 1

    def fix_asl_meta(self, meta, asl):
        if meta.get('RepetitionTimePreparation', 0) > 100:
            meta['RepetitionTimePreparation'] /= 1000
            return ['RepetitionTimePreparation ms -> s']
        return []

    def fix_context(self, context, asl, data):
        n_c, n_l = context.count('control'), context.count('label')
        if n_c == n_l:
            return context, data, []
        extra = 'control' if n_c > n_l else 'label'
        drop = max(i for i, t in enumerate(context) if t == extra)
        keep = [i for i in range(len(context)) if i != drop]
        return [context[i] for i in keep], data[..., keep], [f'dropped unpaired {extra} {drop}']


class PennSpiral(Spec):
    min_sessions = 1

    def fix_asl_meta(self, meta, asl):
        if isinstance(meta.get('TotalAcquiredPairs'), float):
            meta['TotalAcquiredPairs'] = int(meta['TotalAcquiredPairs'])
            return ['TotalAcquiredPairs -> int']
        return []


class DS007116(PennSpiral):
    name, exemplar = 'ds007116', 'sub-20721'
    t1w_pattern = r'_rec-normrefaced_run-01_T1w\.nii\.gz$'
    min_sessions = 2


class DS008579(PennSpiral):
    name, exemplar = 'ds008579', 'sub-112061'
    asl_pattern = r'asl\.nii\.gz'  # the standard protocol only (no CuBIDS variants)
    t1w_pattern = r'_rec-refaced_T1w\.nii\.gz$'


class DS006131(PennSpiral):
    name, exemplar = 'ds006131', 'sub-60295'
    asl_pattern = r'acq-v2_run-01_asl\.nii\.gz'
    t1w_pattern = r'_rec-norm_T1w\.nii\.gz$'


class DS005752(Spec):
    name, exemplar = 'ds005752', 'sub-ON13986'
    t1w_pattern = r'_acq-MPRAGE_rec-SCIC_T1w\.nii\.gz$'
    separate_m0 = False
    root_aslcontext = True
    max_sessions = 1


SPECS = {
    s.name: s
    for s in (DS008719, DS004856, DS004146, DS006148, DS007116, DS008579, DS006131, DS005752)
}


# ---------------------------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------------------------
def choose_subjects(spec, n):
    subs = sorted(p for p in spec.src.glob('sub-*') if p.is_dir() and p.name != spec.exemplar)
    eligible = [s for s in subs if spec.eligible(s)]
    rng = random.Random(f'{SEED}-{spec.name}')
    return sorted(rng.sample(eligible, min(n, len(eligible))), key=lambda p: p.name), len(eligible)


def build(spec, n, out):
    subjects, n_eligible = choose_subjects(spec, n)
    dest = out / spec.name
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    _write_json(
        dest / 'dataset_description.json',
        {
            'Name': f'{spec.name} (HMC probe extension)',
            'BIDSVersion': '1.10.0',
            'DatasetType': 'raw',
        },
    )
    log = [
        f'# {spec.name}: HMC probe extension',
        '',
        f'{len(subjects)} of {n_eligible} eligible subjects (seed {SEED}), from {spec.src}.',
        *[f'- {note}' for note in spec.notes],
        '',
    ]
    for sub_dir in subjects:
        runs = spec.runs(sub_dir)
        files = []
        for ses_dir, asls in runs.items():
            files.append(spec.t1w(ses_dir))
            for asl in asls:
                files.append(asl)
                if spec.m0(asl):
                    files.append(spec.m0(asl))
        _cmd_datalad(spec.src, 'get', files)
        for ses_dir, asls in runs.items():
            rel_ses = ses_dir.relative_to(spec.src)
            t1w = spec.t1w(ses_dir)
            t1_dst = dest / rel_ses / 'anat' / spec.t1w_name(t1w)
            t1_dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(t1w, t1_dst)
            t1_json = t1w.with_name(t1w.name.replace('.nii.gz', '.json'))
            if t1_json.exists():
                shutil.copyfile(t1_json, t1_dst.with_name(t1_dst.name.replace('.nii.gz', '.json')))
            for asl in asls:
                rel = asl.relative_to(spec.src)
                meta = _read_json(asl.with_name(asl.name.replace('.nii.gz', '.json')))
                changes = spec.fix_asl_meta(meta, asl)
                ctx_src = spec.aslcontext(asl)
                context = pd.read_table(ctx_src, keep_default_na=False)['volume_type'].tolist()
                img = nb.load(asl)
                data = np.asanyarray(img.dataobj)
                context, data, more = spec.fix_context(context, asl, data.astype(np.float32))
                changes += more
                (dest / rel).parent.mkdir(parents=True, exist_ok=True)
                if data.shape[-1] != img.shape[-1]:
                    for key, value in list(meta.items()):
                        if (
                            isinstance(value, list)
                            and len(value) == img.shape[-1]
                            and key != 'SliceTiming'
                        ):
                            meta[key] = [v for i, v in enumerate(value) if i < data.shape[-1]]
                    img = nb.Nifti1Image(data, img.affine, img.header)
                    img.to_filename(dest / rel)
                else:
                    shutil.copyfile(asl, dest / rel)
                _write_json(dest / str(rel).replace('.nii.gz', '.json'), meta)
                pd.DataFrame({'volume_type': context}).to_csv(
                    dest / str(rel).replace('_asl.nii.gz', '_aslcontext.tsv'),
                    sep='\t',
                    index=False,
                )
                m0 = spec.m0(asl)
                if m0:
                    m0_rel = m0.relative_to(spec.src)
                    shutil.copyfile(m0, dest / m0_rel)
                    m0_meta = _read_json(m0.with_name(m0.name.replace('.nii.gz', '.json')))
                    m0_meta['IntendedFor'] = [f'bids::{rel.as_posix()}']
                    _write_json(dest / str(m0_rel).replace('.nii.gz', '.json'), m0_meta)
                log.append(f'- {rel}: ' + ('; '.join(changes) if changes else 'copied unchanged'))
        _cmd_datalad(spec.src, 'drop', [f for f in files if f is not None])
    (dest / 'EXTENDED.md').write_text('\n'.join(log) + '\n')
    names = [s.name for s in subjects]
    held_out = sorted(
        random.Random(f'{SEED}-{spec.name}-holdout').sample(names, min(N_HELD_OUT, len(names)))
    )
    return {'subjects': names, 'held_out': held_out, 'n_eligible': n_eligible}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--n', type=int, default=20)
    args = parser.parse_args(argv)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    split_file = out / 'split.json'
    split = json.loads(split_file.read_text()) if split_file.exists() else {}
    for name, cls in SPECS.items():
        if args.only and name not in args.only:
            continue
        split[name] = build(cls(), args.n, out)
        split_file.write_text(json.dumps(split, indent=2) + '\n')
        print(
            f'{name}: {len(split[name]["subjects"])} subjects of {split[name]["n_eligible"]}',
            flush=True,
        )


if __name__ == '__main__':
    main()
