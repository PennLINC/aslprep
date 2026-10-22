#!/usr/bin/env python3
"""Curate the panel-extension subjects into BIDS datasets outside the zoo.

Plan: ``docs/plans/2026-10-22-joint-multicontrast-coregistration.md``, "Panel extension".
Each function adapts the curator of the same dataset in the zoo's ``code/03_curate_openneuro.py``
from one exemplar subject to the drawn subjects (``extension_subjects.tsv``), with the zoo's own
``Builder`` writing to ``--out`` instead of the zoo. Every change from the source is recorded in
the curated dataset's README, as in the zoo. The zoo itself is never written to.

Usage (in WSL)::

    python curate_extension.py --zoo-code /mnt/c/.../aslprep-zoo/code \\
        --subjects extension_subjects.tsv --out /home/tsalo/coreg-mc/panel-ext [ds ...]
"""

import argparse
import csv
import re
import sys
from pathlib import Path

zu = None  # the zoo's zoo_utils, imported in main()


def _asl_sessions(src, sub):
    """Session directories of ``sub`` (relative to ``src``) that contain ASL."""
    sessions = sorted(p for p in (src / sub).glob('ses-*') if p.is_dir()) or [src / sub]
    return [s.relative_to(src).as_posix() for s in sessions if list((s / 'perf').glob('*_asl.nii*'))]


def header(ds, version, subjects):
    return (
        f'Coregistration-probe panel extension from OpenNeuro {ds} (snapshot {version}).\n'
        f'Subjects kept ({len(subjects)}): {", ".join(subjects)}.\n'
        f'Source: https://openneuro.org/datasets/{ds}/versions/{version}\n'
        'Curated by aslprep scripts/coreg_mc/curate_extension.py, adapting the aslprep-zoo curators.'
    )


def curate_ds003823(subjects):
    b = zu.Builder('ds003823', zu.CACHE / 'ds003823')
    b.top_level(subjects)
    runs = {'': 'acq-rest', 'task-trainingMimickingASL_': 'acq-trainingMimicking'}
    for sub in subjects:
        for base in _asl_sessions(b.src, sub):
            ses = base.split('/')[-1]
            t1w = f'{base}/anat/{sub}_{ses}_acq-wholebrain_T1w.nii.gz'
            if (b.src / t1w).exists():
                b.image(t1w, t1w)
            for src_ent, acq in runs.items():
                src_stem = f'{base}/perf/{sub}_{ses}_{src_ent}'
                if not (b.src / f'{src_stem}asl.nii').exists():
                    continue
                dst_stem = f'{base}/perf/{sub}_{ses}_{acq}_'
                n_vols = zu.n_volumes(b.src / f'{src_stem}asl.nii')

                def fix_asl(meta, n_vols=n_vols):
                    meta['TotalAcquiredPairs'] = n_vols // 2
                    meta.pop('PulseSequenceType', None)
                    meta.pop('SliceTiming', None)

                def fix_m0(meta, dst_stem=dst_stem):
                    meta['IntendedFor'] = [f'bids::{dst_stem}asl.nii.gz']
                    for key in ('PulseSequenceType', 'TotalAcquiredPairs', 'M0Type'):
                        meta.pop(key, None)

                b.image(f'{src_stem}asl.nii', f'{dst_stem}asl.nii.gz', edit=fix_asl)
                b.image(f'{src_stem}m0scan.nii', f'{dst_stem}m0scan.nii.gz', edit=fix_m0)
                zu.copy(b.src / f'{src_stem}aslcontext.tsv', b.out / f'{dst_stem}aslcontext.tsv')
    b.note(
        'perf files renamed: the task entity is not allowed for asl/m0scan. "*_asl" -> '
        '"*_acq-rest_asl" and "*_task-trainingMimickingASL_asl" -> "*_acq-trainingMimicking_asl" '
        '(same for m0scan and aslcontext). Physiological recordings were not downloaded.'
    )
    b.note('asl.json: TotalAcquiredPairs set to half the number of volumes; removed the deprecated '
           'PulseSequenceType and the all-zero SliceTiming (3D readout).')
    b.note('m0scan.json: added the missing IntendedFor; removed fields that only apply to asl.')
    b.note('anat: kept the acq-wholebrain T1w only (acq-lc is a locus coeruleus slab).')
    b.note('Uncompressed .nii images were gzipped. Only sessions with ASL were kept.')
    tsv = b.out / 'participants.tsv'
    lines = tsv.read_text().splitlines()
    fixed = [lines[0]] + ['\t'.join('n/a' if c == 'N/A' else c for c in ln.split('\t')) for ln in lines[1:]]
    if fixed != lines:
        zu.write_text(tsv, '\n'.join(fixed) + '\n')
        b.note('participants.tsv: "N/A" cells rewritten as the BIDS missing value "n/a".')
    b.finish(header('ds003823', '1.3.5', subjects))


def curate_ds004146(subjects):
    b = zu.Builder('ds004146', zu.CACHE / 'ds004146')
    b.top_level(subjects)
    for sub in subjects:
        for base in _asl_sessions(b.src, sub):
            ses = base.split('/')[-1]
            stem = f'{base}/anat/{sub}_{ses}'
            if (b.src / f'{stem}_UNIT1_denoised.nii.gz').exists():
                # Copy with the T1w suffix, so that ASLPrep (which requires a T1w) can run.
                b.image(f'{stem}_UNIT1_denoised.nii.gz', f'{stem}_rec-denoised_T1w.nii.gz')
            if (b.src / f'{stem}_T2w.nii.gz').exists():
                b.image(f'{stem}_T2w.nii.gz', f'{stem}_T2w.nii.gz')
            stem = f'{base}/perf/{sub}_{ses}'
            b.image(f'{stem}_asl.nii.gz', f'{stem}_asl.nii.gz')
            zu.write_aslcontext(
                b.out / f'{stem}_aslcontext.tsv', zu.read_aslcontext(b.src / f'{stem}_aslcontext.tsv')
            )
    b.note(
        'anat: the source has no T1w, only MP2RAGE UNIT1. The denoised UNIT1 was copied with the '
        'T1w suffix ("*_rec-denoised_T1w") so that ASLPrep has a T1w; other UNIT1, MP2RAGE, '
        'FLAIR and T2w-TSE images were not downloaded.'
    )
    b.note(
        'perf: the separate m0scan was not kept. M0Type is "Included" and a separate M0 is also '
        'associated, which ASLPrep rejects (the HMC probe adaptation of this dataset).'
    )
    b.note('aslcontext.tsv rewritten with LF line endings and without trailing blank lines.')
    b.note('fmap not kept: the EPI fieldmaps have no IntendedFor and match the BOLD/DWI geometry.')
    pjson = zu.read_json(b.out / 'participants.json')
    rows = [ln.split('\t') for ln in (b.out / 'participants.tsv').read_text().splitlines()]
    for column in ('ses01_protocol', 'ses02_protocol'):
        if column not in rows[0] or column not in pjson:
            continue
        levels = pjson[column].setdefault('Levels', {})
        for row in rows[1:]:
            value = row[rows[0].index(column)]
            if value != 'n/a' and value not in levels:
                levels[value] = 'Not described in the source participants.json'
                b.note(f'participants.json: added the undefined {column} level "{value}".')
    zu.write_json(b.out / 'participants.json', pjson)
    b.finish(header('ds004146', '1.0.0', subjects))


def curate_ds008020(subjects):
    b = zu.Builder('ds008020', zu.CACHE / 'ds008020')
    b.top_level(subjects)
    for sub in subjects:
        unit1 = next((b.src / sub / 'anat').glob(f'{sub}_UNIT1.nii*'))
        b.image(f'{sub}/anat/{unit1.name}', f'{sub}/anat/{sub}_T1w.nii.gz')
        for src in sorted((b.src / sub / 'perf').glob('*_asl.nii')):
            stem = f'{sub}/perf/{src.name[:-8]}'

            def fix(meta):
                if meta['RepetitionTimePreparation'] > 100:  # milliseconds
                    meta['RepetitionTimePreparation'] = meta['RepetitionTimePreparation'] / 1000
                if 'PulseSequenceType' in meta:
                    meta['PulseSequenceDetails'] = meta.pop('PulseSequenceType')

            b.image(f'{stem}_asl.nii', f'{stem}_asl.nii.gz', edit=fix)
            zu.copy(b.src / f'{stem}_aslcontext.tsv', b.out / f'{stem}_aslcontext.tsv')
    b.note('asl.json: RepetitionTimePreparation converted from milliseconds to seconds where it '
           'was above 100; deprecated PulseSequenceType moved to PulseSequenceDetails.')
    b.note('anat: the source has no T1w, only a brain-extracted MP2RAGE UNIT1, copied with the T1w '
           'suffix ("*_T1w", identical data, still brain-extracted). Run with --skull-strip-t1w skip.')
    b.note('All images written as .nii.gz. Per-subject derivatives/ folders not copied.')
    b.finish(header('ds008020', '1.0.4', subjects))


def curate_ds007522(subjects):
    b = zu.Builder('ds007522', zu.CACHE / 'ds007522')
    b.top_level(subjects)
    for sub in subjects:
        for base in _asl_sessions(b.src, sub):
            anat = sorted(p.name for p in (b.src / base / 'anat').glob('*.nii.gz'))
            has_t2w = any(n.endswith('_T2w.nii.gz') for n in anat)
            for name in anat:
                if name.endswith(('_T1w.nii.gz', '_T2w.nii.gz')):
                    b.image(f'{base}/anat/{name}', f'{base}/anat/{name}')
                elif name.endswith('_FLAIR.nii.gz') and not has_t2w:
                    b.image(f'{base}/anat/{name}', f'{base}/anat/{name}')
            perf = sorted(p.name for p in (b.src / base / 'perf').glob('*.nii.gz'))
            asl_names = [n for n in perf if n.endswith('_asl.nii.gz')]
            ses = base.split('/')[-1]
            for name in perf:
                if 'RR_' in name:
                    continue  # duplicates of the acq-pcasl3000 M0 scans

                def fix(meta, name=name, ses=ses, asl_names=asl_names):
                    if meta.get('ArterialSpinLabelingType') == 'PASL':
                        meta.pop('LabelingDuration', None)
                        meta.pop('BolusCutOffDelayTechnique', None)
                        meta.pop('PostLabelDelay', None)
                    if name.endswith('_m0scan.nii.gz') and '_dir-' in name:
                        direction = re.search(r'_(dir-[A-Z]+)_', name).group(1)
                        meta['IntendedFor'] = [
                            f'{ses}/perf/{n}' for n in asl_names if f'_{direction}_' in n
                        ]

                b.image(f'{base}/perf/{name}', f'{base}/perf/{name}', edit=fix)
                if name.endswith('_asl.nii.gz'):
                    ctx = name.replace('_asl.nii.gz', '_aslcontext.tsv')
                    zu.copy(b.src / base / 'perf' / ctx, b.out / base / 'perf' / ctx)
            for path in sorted((b.src / base / 'fmap').glob('*.nii.gz')):
                b.image(f'{base}/fmap/{path.name}', f'{base}/fmap/{path.name}')
    b.prune_intended_for()
    b.note('acq-pcasl3000RR m0scans dropped (the same series as the acq-pcasl3000 m0scans).')
    b.note('m0scan IntendedFor narrowed to the runs with the same phase-encoding direction.')
    b.note('PASL asl.json: removed LabelingDuration, BolusCutOffDelayTechnique and PostLabelDelay.')
    b.note('Fieldmap IntendedFor lists pruned to the files kept; fieldmaps left without targets dropped.')
    b.note('anat: T1w plus T2w, or 3D FLAIR where there is no T2w. Only sessions with ASL kept.')
    b.finish(header('ds007522', '2.0.1', subjects))


def curate_ds007116(subjects):
    b = zu.Builder('ds007116', zu.CACHE / 'ds007116')
    b.top_level(subjects)
    for sub in subjects:
        for base in _asl_sessions(b.src, sub):
            for path in sorted((b.src / base / 'anat').glob('*.nii.gz')):
                is_t1w = '_rec-normrefaced_' in path.name and path.name.endswith('_T1w.nii.gz')
                is_t2w = '_rec-normdefaced_' in path.name and path.name.endswith('_T2w.nii.gz')
                if is_t1w or is_t2w:
                    b.image(f'{base}/anat/{path.name}', f'{base}/anat/{path.name}')

            def fix(meta):
                pairs = meta.get('TotalAcquiredPairs')
                if isinstance(pairs, float):
                    meta['TotalAcquiredPairs'] = int(pairs)

            for path in sorted((b.src / base / 'perf').glob('*')):
                rel = f'{base}/perf/{path.name}'
                if path.name.endswith('.nii.gz'):
                    b.image(rel, rel, edit=fix)
                elif path.name.endswith('_aslcontext.tsv'):
                    zu.copy(path, b.out / rel)
    b.prune_intended_for()
    b.note('anat: kept the rec-normrefaced T1w and rec-normdefaced T2w.')
    b.note('Fieldmaps not kept: they are IntendedFor BOLD or DWI runs.')
    b.note('perf file names (including CuBIDS "acq-VARIANT..." labels) and aslcontext kept as in '
           'the source; TotalAcquiredPairs made an integer where it was a float.')
    b.finish(header('ds007116', '1.0.6', subjects))


def curate_ds005752(subjects):
    b = zu.Builder('ds005752', zu.CACHE / 'ds005752')
    b.top_level(subjects)
    context = zu.read_aslcontext(b.src / 'aslcontext.tsv')
    for sub in subjects:
        for base in _asl_sessions(b.src, sub):
            ses = base.split('/')[-1]
            anat = sorted((b.src / base / 'anat').glob('*.nii.gz'))
            for suffix in ('T1w', 'T2w'):
                files = [p for p in anat if p.name.endswith(f'_{suffix}.nii.gz')]
                scic = [p for p in files if '_rec-SCIC_' in p.name]
                for path in scic or files:
                    b.image(f'{base}/anat/{path.name}', f'{base}/anat/{path.name}')
            stem = f'{base}/perf/{sub}_{ses}'
            asl = b.image(f'{stem}_asl.nii.gz', f'{stem}_asl.nii.gz')
            assert zu.n_volumes(asl) == len(context), asl
            zu.write_aslcontext(b.out / f'{stem}_aslcontext.tsv', context)
    b.note('aslcontext.tsv (m0scan, deltam, cbf) copied from the dataset root to sit next to each '
           'asl file, as in the zoo exemplar.')
    b.note('anat: one T1w and one T2w per session, preferring rec-SCIC; FLAIR not kept.')
    b.finish(header('ds005752', '2.1.0', subjects))


CURATORS = {
    'ds007116': curate_ds007116,
    'ds003823': curate_ds003823,
    'ds004146': curate_ds004146,
    'ds008020': curate_ds008020,
    'ds007522': curate_ds007522,
    'ds005752': curate_ds005752,
}


def main():
    global zu
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--zoo-code', type=Path, required=True)
    parser.add_argument('--subjects', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('datasets', nargs='*')
    args = parser.parse_args()

    out = args.out.resolve()
    zoo = args.zoo_code.resolve().parent
    # Builder rebuilds its output directory from scratch, so the root must be a dedicated one.
    if not out.is_absolute() or out == out.anchor or zoo in (out, *out.parents) or out.name != 'panel-ext':
        sys.exit(f'refusing output root {out}')
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.zoo_code))
    sys.dont_write_bytecode = True  # importing must not write __pycache__ into the zoo
    import zoo_utils

    zoo_utils.ZOO = out  # Builder writes to ZOO / <dataset>
    zu = zoo_utils

    with args.subjects.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    for ds in args.datasets or CURATORS:
        subjects = [r['subject'] for r in rows if r['dataset'] == ds]
        CURATORS[ds](subjects)


if __name__ == '__main__':
    main()
