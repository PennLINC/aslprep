"""Shared machinery for the truth-scored integration tests (spec Section 4.8).

Each ``test_aslscan_<recipe>.py`` module runs ASLPrep once on its fixture (a module-scoped
fixture built with :func:`run_recipe`), writes ``truth_score.json`` before any assertion, and
then checks one metric per test item, so one failure does not hide the others.
:func:`shared_items` adds the invariants every module must check.
"""

import functools
import json
import math
import shutil
from dataclasses import dataclass
from pathlib import Path

import pytest

from aslprep.tests import aslscan_fixtures as af
from aslprep.tests import truth_bounds as tb
from aslprep.tests import truth_scoring as ts


@dataclass(frozen=True)
class Known:
    """A reported, expected failure, held to the size it was diagnosed at.

    The item must fail, and the metric at ``path`` must lie in ``[lo, hi]``: a fix (the item
    passes) or a change in size (worse or better) fails the item instead of hiding behind the
    expected failure.
    """

    reason: str
    path: tuple
    lo: float
    hi: float


def known_failure(fixture_name, known):
    """Decorate a test item taking ``request`` as a :class:`Known` failure."""

    def decorate(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            try:
                fn(*args, **kwargs)
            except AssertionError as exc:
                score = kwargs['request'].getfixturevalue(fixture_name).score
                value = ts.lookup(score, known.path)
                label = '.'.join(known.path)
                if not _finite(value) or not known.lo <= value <= known.hi:
                    pytest.fail(
                        f'known failure changed size: {label} = {value!r}, known range '
                        f'[{known.lo}, {known.hi}] ({known.reason}); the check said: {exc}'
                    )
                pytest.xfail(f'{known.reason} [{label} = {value:.4g}]')
            pytest.fail(f'known failure now passes; remove it: {known.reason}')

        return wrapper

    return decorate


def _finite(value):
    return isinstance(value, int | float) and math.isfinite(value)


@dataclass
class RecipeRun:
    recipe: str
    test_name: str
    fixture: Path
    out_dir: Path
    score: dict


def run_recipe(
    recipe_name,
    test_name,
    data_dir,
    output_dir,
    working_dir,
    spaces=('asl',),
    extra_args=(),
    fwhm=5.0,
    score_spaces=(),
    extra_cbf=(),
):
    """Run ASLPrep on a fixture, score it against the truth, and write ``truth_score.json``.

    A pipeline failure is recorded in the score (``error``) and re-raised after writing it, so
    the artifact always exists.
    """
    from aslprep.tests.test_cli import _run_and_generate

    recipe = af.load_recipe(recipe_name)
    fixture = af.fixture_dir(recipe_name, data_dir)
    test_dir = Path(output_dir) / test_name
    out_dir = test_dir / 'aslprep'
    work_dir = Path(working_dir) / test_name
    # Start clean: files left by an earlier local run would leak into the output manifest.
    shutil.rmtree(test_dir, ignore_errors=True)
    shutil.rmtree(work_dir, ignore_errors=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    filter_file = test_dir / 'bids_filters.json'
    filter_file.write_text(json.dumps({'asl': {'acquisition': recipe.acq}}, indent=2))

    parameters = [
        str(fixture),
        str(out_dir),
        'participant',
        '--participant-label=01',
        f'-w={work_dir}',
        f'--bids-filter-file={filter_file}',
        '--output-spaces',
        *spaces,
        '--fs-no-reconall',
        f'--smooth_kernel={fwhm}',
    ]
    if recipe.anat == 'derivatives':
        parameters += ['--derivatives', f'anat={fixture / "derivatives" / "anat"}']
    if recipe.m0_divisor != 1:
        parameters.append(f'--m0_scale={recipe.m0_divisor}')
    parameters += list(extra_args)

    try:
        _run_and_generate(test_name, '01', parameters, str(out_dir), check_outputs=False)
    except Exception as exc:  # recorded in the score, then re-raised
        _write_score(test_dir, {'recipe': recipe_name, 'error': f'{type(exc).__name__}: {exc}'})
        raise
    score = ts.score_run(
        fixture, out_dir, recipe.acq, fwhm=fwhm, spaces=score_spaces, extra_cbf=extra_cbf
    )
    _write_score(test_dir, score)
    return RecipeRun(recipe_name, test_name, fixture, out_dir, score)


def _write_score(test_dir, score):
    path = Path(test_dir) / 'truth_score.json'
    path.write_text(json.dumps(score, indent=2, sort_keys=True, default=_jsonable) + '\n')
    print(f'TRUTH SCORE {Path(test_dir).name} {json.dumps(score, default=_jsonable)}')


def _jsonable(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return str(value)


def shared_items(fixture_name, multi_delay=False, report_only=(), known=None):
    """Test functions every integration module includes (spec Section 4.7, metric contract).

    ``fixture_name`` names the module's run fixture. ``report_only`` lists score paths that
    are not asserted against bounds but must exist and be finite. ``known`` maps item names
    to a :class:`Known` failure.
    """

    def _run(request):
        return request.getfixturevalue(fixture_name)

    def test_clean_run(request):
        clean = _run(request).score['clean']
        if clean['crash_files']:
            pytest.fail(f'crash files: {clean["crash_files"]}')
        if not clean['report_ok']:
            pytest.fail('the HTML report lists errors')

    def test_outputs_manifest(request):
        from aslprep.tests.utils import check_generated_files, get_test_data_path

        run = _run(request)
        manifest = Path(get_test_data_path()) / f'expected_outputs_{run.test_name}.txt'
        if not manifest.exists():
            found = sorted(
                p.relative_to(run.out_dir).as_posix()
                for p in run.out_dir.rglob('*')
                if p.is_file() and 'figures' not in p.parts and 'log' not in p.parts
            )
            candidate = run.out_dir.parent / manifest.name
            candidate.write_text('\n'.join(found) + '\n')
            pytest.fail(f'{manifest.name} does not exist; a candidate was written to {candidate}')
        check_generated_files(str(run.out_dir), str(manifest))

    def test_metric_contract(request):
        run = _run(request)
        missing = [
            '.'.join(path) for path in report_only if not _finite(ts.lookup(run.score, path))
        ]
        problems = [
            p
            for section in run.score.values()
            if isinstance(section, dict)
            for p in section.get('problems', [])
        ]
        if missing:
            pytest.fail(f'report-only metrics missing or not finite: {missing}')
        if problems:
            pytest.fail(f'scoring problems: {problems}')

    def test_mask_coverage(request):
        run = _run(request)
        tb.check(run.score, 'mask_coverage', run.recipe)

    def test_finite(request):
        run = _run(request)
        tb.check(run.score, 'finite_multi' if multi_delay else 'finite', run.recipe)

    @pytest.mark.parametrize('tissue', ['GM', 'WM'])
    def test_n_dominant(request, tissue):
        run = _run(request)
        tb.check(run.score, 'n_dominant', run.recipe, t=tissue)

    items = {
        'test_clean_run': test_clean_run,
        'test_outputs_manifest': test_outputs_manifest,
        'test_metric_contract': test_metric_contract,
        'test_mask_coverage': test_mask_coverage,
        'test_finite': test_finite,
        'test_n_dominant': test_n_dominant,
    }

    if multi_delay:

        def test_fit_bound(request):
            """Brain voxels whose multi-delay fit ends at the CBF bound have no estimate."""
            run = _run(request)
            tb.check(run.score, 'fit_bound', run.recipe)

        items['test_fit_bound'] = test_fit_bound
    for name, entry in (known or {}).items():
        items[name] = known_failure(fixture_name, entry)(items[name])
    return items


def scored_items(fixture_name, recipe, spaces=(), suppression_pulses=None, known=None):
    """The truth-scored items common to the no-motion recipes (spec Sections 4.7 and 5).

    ``spaces`` are score keys such as ``'space-MNI152NLin2009cAsym'``. ``suppression_pulses``
    adds the physical-agreement item for background suppression (spec Section 8). ``known``
    maps item names to a :class:`Known` failure.
    """

    def _score(request):
        return request.getfixturevalue(fixture_name).score

    def test_tier_a_median(request):
        tb.check(_score(request), 'tier_a_median', recipe)

    def test_tier_a_quantification(request):
        run = request.getfixturevalue(fixture_name)
        score = run.score
        tb.check(score, 'tier_a_quant_median', recipe)
        if ts.is_multi_delay(_asl_sidecar(run)):
            tb.check(score, 'tier_a_quant_agree', recipe)
            tb.check(score, 'tier_a_att', recipe)
        else:
            tb.check(score, 'tier_a_quant_p95', recipe)

    def test_motion_correction(request):
        """Without motion, motion correction must keep every volume in register with the rest."""
        score = _score(request)
        tb.check(score, 'motion_median', recipe)
        tb.check(score, 'motion_max', recipe)

    @pytest.mark.parametrize('tissue', ['GM', 'WM'])
    def test_tier_b(request, tissue):
        score = _score(request)
        reference = ('native', 'expected_ratio', '{t}')
        # multi-delay expectations exist on a voxel subset: compare on the same voxels
        name = 'tier_b_matched' if ts.lookup(score, ('native', 'tier_b_matched')) else 'tier_b'
        tb.check(score, name, recipe, reference=reference, t=tissue)

    def test_coregistration(request):
        score = _score(request)
        tb.check(score, 'coreg_rot', recipe)
        tb.check(score, 'coreg_rms', recipe)

    items = {
        'test_tier_a_median': test_tier_a_median,
        'test_tier_a_quantification': test_tier_a_quantification,
        'test_motion_correction': test_motion_correction,
        'test_tier_b': test_tier_b,
        'test_coregistration': test_coregistration,
    }

    if spaces:

        @pytest.mark.parametrize('space', list(spaces))
        def test_standard_space(request, space):
            score = _score(request)
            tb.check(score, 'space_alignment', recipe, space=space)
            reference = ('native', 'expected_ratio', '{t}')
            tb.check(score, 'space', recipe, reference=reference, space=space, t='GM')
            if space == 'space-T1w':
                tb.check(score, 'space_resampling', recipe, space=space)
            for tissue in ('GM', 'WM'):
                tb.check(score, 'space_finite', recipe, space=space, t=tissue)

        items['test_standard_space'] = test_standard_space

    if suppression_pulses:

        def test_background_suppression(request):
            """ASLPrep's suppression loss against the simulator's (spec Section 8).

            Both keep about 95 % of the ASL signal per pulse (Alsop et al. 2015). Without
            LabelingEfficiency in the sidecar, ASLPrep applies 0.95^n and must match the
            physically calibrated CBF. With it, ASLPrep applies no suppression loss, so its CBF
            is 0.95^n of the physical value (#706): asserted at that size, then marked expected.
            """
            run = request.getfixturevalue(fixture_name)
            score = run.score
            if 'LabelingEfficiency' not in _asl_sidecar(run):
                tb.check(score, 'physical_median', recipe)
                return
            predicted = 0.95**suppression_pulses
            measured = ts.lookup(score, ('native', 'physical', 'median'))
            if measured is None or abs(measured / predicted - 1) > 0.05:
                pytest.fail(
                    f'the suppression result changed: measured {measured}, '
                    f'predicted {predicted:.4f}'
                )
            pytest.xfail(
                'PennLINC/aslprep#706: with LabelingEfficiency in the sidecar ASLPrep applies no '
                f'suppression loss ({suppression_pulses} pulses), so CBF is {measured:.4f} of the '
                f'physically calibrated value (predicted {predicted:.4f})'
            )

        items['test_background_suppression'] = test_background_suppression
    for name, entry in (known or {}).items():
        items[name] = known_failure(fixture_name, entry)(items[name])
    return items


def _asl_sidecar(run):
    acq = af.load_recipe(run.recipe).acq
    return json.loads((run.fixture / 'sub-01' / 'perf' / f'sub-01_acq-{acq}_asl.json').read_text())
