# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Parser."""

import sys

from aslprep import config


def _build_parser():
    """Build parser object."""
    from argparse import Action, ArgumentDefaultsHelpFormatter, ArgumentParser
    from functools import partial
    from pathlib import Path

    from niworkflows.utils.spaces import OutputReferencesAction, Reference
    from packaging.version import Version

    from aslprep.cli.version import check_latest, is_flagged

    deprecations = {
        # parser attribute name: (replacement flag, version slated to be removed in)
        'force_bbr': ('--force bbr', '26.0.0'),
        'force_no_bbr': ('--force no-bbr', '26.0.0'),
        'force_syn': ('--force syn-sdc', '26.0.0'),
    }

    class DeprecatedAction(Action):
        def __call__(self, parser, namespace, values, option_string=None):
            new_opt, rem_vers = deprecations.get(self.dest, (None, None))
            msg = (
                f'{self.option_strings} has been deprecated and will be removed in '
                f'{rem_vers or "a later version"}.'
            )
            if new_opt:
                msg += f' Please use `{new_opt}` instead.'
            print(msg, file=sys.stderr)
            delattr(namespace, self.dest)

    class ToDict(Action):
        def __call__(self, parser, namespace, values, option_string=None):
            d = {}
            for spec in values:
                try:
                    name, loc = spec.split('=')
                    loc = Path(loc)
                except ValueError:
                    loc = Path(spec)
                    name = loc.name

                if name in d:
                    raise ValueError(f'Received duplicate derivative name: {name}')

                d[name] = loc
            setattr(namespace, self.dest, d)

    def _path_exists(path, parser):
        """Ensure a given path exists."""
        if path is None or not Path(path).exists():
            raise parser.error(f'Path does not exist: <{path}>.')
        return Path(path).absolute()

    def _is_file(path, parser):
        """Ensure a given path exists and it is a file."""
        path = _path_exists(path, parser)
        if not path.is_file():
            raise parser.error(f'Path should point to a file (or symlink of file): <{path}>.')
        return path

    def _min_one(value, parser):
        """Ensure an argument is not lower than 1."""
        value = int(value)
        if value < 1:
            raise parser.error("Argument can't be less than one.")
        return value

    def _to_gb(value):
        scale = {'G': 1, 'T': 10**3, 'M': 1e-3, 'K': 1e-6, 'B': 1e-9}
        digits = ''.join([c for c in value if c.isdigit()])
        units = value[len(digits) :] or 'M'
        return int(digits) * scale[units[0]]

    def _process_value(value):
        import bids

        if value is None:
            return bids.layout.Query.NONE
        elif value == '*':
            return bids.layout.Query.ANY
        else:
            return value

    def _filter_pybids_none_any(dct):
        d = {}
        for k, v in dct.items():
            if isinstance(v, list):
                d[k] = [_process_value(val) for val in v]
            else:
                d[k] = _process_value(v)
        return d

    def _bids_filter(value, parser):
        from json import JSONDecodeError, loads

        if value:
            if Path(value).exists():
                try:
                    return loads(Path(value).read_text(), object_hook=_filter_pybids_none_any)
                except JSONDecodeError as e:
                    raise parser.error(f'JSON syntax error in: <{value}>.') from e
            else:
                raise parser.error(f'Path does not exist: <{value}>.')

    verstr = f'ASLPrep v{config.environment.version}'
    currentv = Version(config.environment.version)
    is_release = not any((currentv.is_devrelease, currentv.is_prerelease, currentv.is_postrelease))

    parser = ArgumentParser(
        description=f'ASLPrep: ASL PREProcessing workflows v{config.environment.version}',
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    PathExists = partial(_path_exists, parser=parser)
    IsFile = partial(_is_file, parser=parser)
    PositiveInt = partial(_min_one, parser=parser)
    BIDSFilter = partial(_bids_filter, parser=parser)

    # Arguments as specified by BIDS-Apps
    # required, positional arguments
    # IMPORTANT: they must go directly with the parser object
    parser.add_argument(
        'bids_dir',
        action='store',
        type=PathExists,
        help=(
            'the root folder of a BIDS valid dataset (sub-XXXXX folders should '
            'be found at the top level in this folder).'
        ),
    )
    parser.add_argument(
        'output_dir',
        action='store',
        type=Path,
        help='the output path for the outcomes of preprocessing and visual reports',
    )
    parser.add_argument(
        'analysis_level',
        choices=['participant'],
        help=(
            'processing stage to be run, only "participant" in the case of '
            'ASLPREP (see BIDS-Apps specification).'
        ),
    )

    # optional arguments
    g_bids = parser.add_argument_group('Options for filtering BIDS queries')
    g_bids.add_argument(
        '--skip-bids-validation',
        '--skip_bids_validation',
        action='store_true',
        default=False,
        help='assume the input dataset is BIDS compliant and skip the validation',
    )
    g_bids.add_argument(
        '--participant-label',
        '--participant_label',
        action='store',
        nargs='+',
        type=lambda label: label.removeprefix('sub-'),
        help=(
            'A space delimited list of participant identifiers or a single identifier '
            '(the sub- prefix can be removed)'
        ),
    )
    # Re-enable when option is actually implemented
    # g_bids.add_argument('-s', '--session-id', action='store', default='single_session',
    #                     help='select a specific session to be processed')
    # Re-enable when option is actually implemented
    # g_bids.add_argument('-r', '--run-id', action='store', default='single_run',
    #                     help='select a specific run to be processed')
    # g_bids.add_argument(
    # "-t", "--task-id", action="store", help="select a specific task to be processed"
    # )
    g_bids.add_argument(
        '--bids-filter-file',
        dest='bids_filters',
        action='store',
        type=BIDSFilter,
        metavar='FILE',
        help=(
            'A JSON file describing custom BIDS input filters using PyBIDS. '
            'For further details, please check out '
            'https://aslprep.readthedocs.io/en/'
            f'{currentv.base_version if is_release else "latest"}/faq.html#'
            'how-do-I-select-only-certain-files-to-be-input-to-ASLPrep'
        ),
    )
    g_bids.add_argument(
        '-d',
        '--derivatives',
        action=ToDict,
        metavar='PACKAGE=PATH',
        type=str,
        nargs='+',
        help=(
            'Search PATH(s) for pre-computed BIDS derivatives. '
            'These may be provided as named folders '
            '(e.g., `--derivatives smriprep=/path/to/smriprep`). '
            'Note: This does not apply to non-BIDS data, like FreeSurfer outputs. '
            'To reuse FreeSurfer outputs, use `--fs-subjects-dir` instead.'
        ),
    )
    g_bids.add_argument(
        '--bids-database-dir',
        metavar='PATH',
        type=Path,
        help=(
            'Path to a PyBIDS database folder, for faster indexing '
            '(especially useful for large datasets). '
            'Will be created if not present.'
        ),
    )

    g_perfm = parser.add_argument_group('Options to handle performance')
    g_perfm.add_argument(
        '--nprocs',
        '--nthreads',
        '--n_cpus',
        '--n-cpus',
        dest='nprocs',
        action='store',
        type=PositiveInt,
        help='maximum number of threads across all processes',
    )
    g_perfm.add_argument(
        '--omp-nthreads',
        action='store',
        type=PositiveInt,
        help='maximum number of threads per-process',
    )
    g_perfm.add_argument(
        '--mem',
        '--mem_mb',
        '--mem-mb',
        dest='memory_gb',
        action='store',
        type=_to_gb,
        help='upper bound memory limit for ASLPrep processes',
    )
    g_perfm.add_argument(
        '--low-mem',
        action='store_true',
        help='attempt to reduce memory usage (will increase disk usage in working directory)',
    )
    g_perfm.add_argument(
        '--use-plugin',
        '--nipype-plugin-file',
        action='store',
        metavar='FILE',
        type=IsFile,
        help='nipype plugin configuration file',
    )
    g_perfm.add_argument(
        '--sloppy',
        action='store_true',
        default=False,
        help='Use low-quality tools for speed - TESTING ONLY',
    )

    g_subset = parser.add_argument_group('Options for performing only a subset of the workflow')
    g_subset.add_argument('--anat-only', action='store_true', help='run anatomical workflows only')
    g_subset.add_argument(
        '--level',
        action='store',
        default='full',
        choices=['minimal', 'resampling', 'full'],
        help=(
            "Processing level; may be 'minimal' (nothing that can be recomputed), "
            "'resampling' (recomputable targets that aid in resampling) "
            "or 'full' (all target outputs). "
            "Currently, 'minimal' and 'resampling' are effectively the same."
        ),
    )
    g_subset.add_argument(
        '--boilerplate-only',
        '--boilerplate_only',
        action='store_true',
        default=False,
        help='generate boilerplate only',
    )
    g_subset.add_argument(
        '--reports-only',
        action='store_true',
        default=False,
        help=(
            "only generate reports, don't run workflows. This will only rerun report "
            'aggregation, not reportlet generation for specific nodes.'
        ),
    )

    g_conf = parser.add_argument_group('Workflow configuration')
    g_conf.add_argument(
        '--ignore',
        required=False,
        action='store',
        nargs='+',
        default=[],
        choices=['fieldmaps', 'sbref', 't2w', 'flair', 'fmap-jacobian'],
        help=(
            'ignore selected aspects of the input dataset to disable corresponding '
            'parts of the workflow (a space delimited list)'
        ),
    )
    g_conf.add_argument(
        '--force',
        required=False,
        action='store',
        nargs='+',
        default=[],
        choices=['bbr', 'no-bbr', 'syn-sdc', 'fmap-jacobian', 'ge', 'no-ge'],
        help='Force selected processing choices, overriding automatic selections '
        '(a space delimited list).\n'
        ' * [no-]bbr: Use/disable boundary-based registration for BOLD-to-T1w coregistration\n'
        '             (No goodness-of-fit checks)\n'
        ' * syn-sdc: Calculate SyN-SDC correction *in addition* to other fieldmaps\n'
        ' * [no-]ge: Use/disable GE-specific processing\n',
    )
    g_conf.add_argument(
        '--output-spaces',
        nargs='*',
        action=OutputReferencesAction,
        help="""\
Standard and non-standard spaces to resample anatomical and functional images to. \
Standard spaces may be specified by the form \
``<SPACE>[:cohort-<label>][:res-<resolution>][...]``, where ``<SPACE>`` is \
a keyword designating a spatial reference, and may be followed by optional, \
colon-separated parameters. \
Non-standard spaces imply specific orientations and sampling grids. \
Important to note, the ``res-*`` modifier does not define the resolution used for \
the spatial normalization. To generate no ASL outputs, use this option without specifying \
any spatial references.""",
    )
    g_conf.add_argument(
        '--longitudinal',
        action='store_true',
        help='treat dataset as longitudinal - may increase runtime',
    )

    g_conf.add_argument(
        '--asl2anat-init',
        action='store',
        choices=['auto', 't1w', 't2w', 'header'],
        default='auto',
        help=(
            'Method of initial ASL to anatomical coregistration. If `auto`, a T2w image is used '
            'if available, otherwise the T1w image. `t1w` forces use of the T1w, `t2w` forces use '
            'of the T2w, and `header` uses the BOLD header information without an initial '
            'registration.'
        ),
    )
    g_conf.add_argument(
        '--asl2anat-dof',
        action='store',
        default=6,
        choices=[6, 9, 12],
        type=int,
        help=(
            'Degrees of freedom when registering ASL to anatomical images. '
            '6 degrees (rotation and translation) are used by default.'
        ),
    )

    g_use_bbr = g_conf.add_mutually_exclusive_group()
    g_use_bbr.add_argument(
        '--force-bbr',
        action=DeprecatedAction,
        help='Deprecated - use `--force bbr` instead. See `--force` for more details.',
    )
    g_use_bbr.add_argument(
        '--force-no-bbr',
        action=DeprecatedAction,
        help='Deprecated - use `--force no-bbr` instead. See `--force` for more details.',
    )
    g_conf.add_argument(
        '--random-seed',
        dest='_random_seed',
        action='store',
        type=int,
        default=None,
        help='Initialize the random seed for the workflow',
    )

    g_use_ge = g_conf.add_mutually_exclusive_group()
    g_use_ge.add_argument(
        '--force-ge',
        action=DeprecatedAction,
        help='Deprecated - use `--force ge` instead. See `--force` for more details.',
    )
    g_use_ge.add_argument(
        '--force-no-ge',
        action=DeprecatedAction,
        help='Deprecated - use `--force no-ge` instead. See `--force` for more details.',
    )

    g_conf.add_argument(
        '--m0-scale',
        '--m0_scale',
        required=False,
        action='store',
        default=1,
        type=float,
        help=(
            'Relative scale between ASL and M0. '
            'M0 scans are multiplied by m0_scale before calculating CBF. '
            'It is important to note, however, that BIDS expects ASL and M0 data to scaled in '
            'the raw dataset, so this parameter should only be used if your dataset does not '
            'have pre-scaled data.'
        ),
    )
    g_conf.add_argument(
        '--smooth-kernel',
        '--smooth_kernel',
        action='store',
        default=5,
        type=float,
        help='Smoothing kernel for the M0 image(s)',
    )
    g_conf.add_argument(
        '--scorescrub',
        action='store_true',
        default=False,
        help="Apply Sudipto Dolui's algorithms for denoising CBF",
    )
    g_conf.add_argument(
        '--basil',
        action='store_true',
        default=False,
        help="FSL's CBF computation with spatial regularization and partial volume correction",
    )

    g_outputs = parser.add_argument_group('Options for modulating outputs')
    g_outputs.add_argument(
        '--aggregate-session-reports',
        dest='aggr_ses_reports',
        action='store',
        type=PositiveInt,
        default=4,
        help="Maximum number of sessions aggregated in one subject's visual report. "
        'If exceeded, visual reports are split by session.',
    )
    g_outputs.add_argument(
        '--medial-surface-nan',
        required=False,
        action='store_true',
        default=False,
        help=(
            'Replace medial wall values with NaNs on functional GIFTI files. '
            'Only performed for GIFTI files mapped to a freesurfer subject '
            '(fsaverage or fsnative).'
        ),
    )
    g_conf.add_argument(
        '--project-goodvoxels',
        required=False,
        action='store_true',
        default=False,
        help=(
            'Exclude voxels whose timeseries have locally high coefficient of variation '
            'from surface resampling. '
            'Only performed for GIFTI files mapped to a freesurfer subject '
            '(fsaverage or fsnative).'
        ),
    )
    g_outputs.add_argument(
        '--md-only-boilerplate',
        action='store_true',
        default=False,
        help='Skip generation of HTML and LaTeX formatted citation with pandoc',
    )
    g_outputs.add_argument(
        '--cifti-output',
        nargs='?',
        const='91k',
        default=False,
        choices=('91k', '170k'),
        type=str,
        help=(
            'Output preprocessed BOLD as a CIFTI dense timeseries. '
            'Optionally, the number of grayordinate can be specified '
            '(default is 91k, which equates to 2mm resolution)'
        ),
    )
    g_outputs.add_argument(
        '--no-msm',
        action='store_false',
        dest='run_msmsulc',
        help='Disable Multimodal Surface Matching surface registration.',
    )

    #  ANTs options
    g_ants = parser.add_argument_group('Specific options for ANTs registrations')
    g_ants.add_argument(
        '--skull-strip-template',
        default='OASIS30ANTs',
        type=Reference.from_string,
        help='select a template for skull-stripping with antsBrainExtraction',
    )
    g_ants.add_argument(
        '--skull-strip-fixed-seed',
        action='store_true',
        help=(
            'do not use a random seed for skull-stripping - will ensure '
            'run-to-run replicability when used with --omp-nthreads 1 and '
            'matching --random-seed <int>'
        ),
    )
    g_ants.add_argument(
        '--skull-strip-t1w',
        action='store',
        choices=('auto', 'skip', 'force'),
        default='force',
        help=(
            "determiner for T1-weighted skull stripping ('force' ensures skull "
            "stripping, 'skip' ignores skull stripping, and 'auto' applies brain extraction "
            'based on the outcome of a heuristic to check whether the brain is already masked).'
        ),
    )

    # Fieldmap options
    g_fmap = parser.add_argument_group('Specific options for handling fieldmaps')
    g_fmap.add_argument(
        '--fmap-bspline',
        action='store_true',
        default=False,
        help='fit a B-Spline field using least-squares (experimental)',
    )
    g_fmap.add_argument(
        '--fmap-no-demean',
        action='store_false',
        default=True,
        help='do not remove median (within mask) from fieldmap',
    )

    # SyN-unwarp options
    g_syn = parser.add_argument_group('Specific options for SyN distortion correction')
    g_syn.add_argument(
        '--use-syn-sdc',
        action='store_true',
        default=False,
        help='EXPERIMENTAL: Use fieldmap-free distortion correction',
    )
    g_syn.add_argument(
        '--force-syn',
        action=DeprecatedAction,
        default=False,
        help='Deprecated - use `--force syn-sdc` instead. See `--force` for more details.',
    )

    # FreeSurfer options
    g_fs = parser.add_argument_group('Specific options for FreeSurfer preprocessing')
    g_fs.add_argument(
        '--fs-license-file',
        metavar='FILE',
        type=IsFile,
        help=(
            'Path to FreeSurfer license key file. Get it (for free) by registering '
            'at https://surfer.nmr.mgh.harvard.edu/registration.html'
        ),
    )
    g_fs.add_argument(
        '--fs-subjects-dir',
        metavar='PATH',
        type=Path,
        help=(
            'Path to existing FreeSurfer subjects directory to reuse. '
            '(default: OUTPUT_DIR/freesurfer)'
        ),
    )
    g_fs.add_argument(
        '--no-submm-recon',
        action='store_false',
        dest='hires',
        help='Disable sub-millimeter (hires) reconstruction',
    )
    g_fs.add_argument(
        '--fs-no-reconall',
        action='store_false',
        dest='run_reconall',
        help='Disable FreeSurfer surface preprocessing.',
    )
    g_fs.add_argument(
        '--fs-no-resume',
        action='store_true',
        dest='fs_no_resume',
        help='EXPERT: Import pre-computed FreeSurfer reconstruction without resuming. '
        'The user is responsible for ensuring that all necessary files are present.',
    )

    g_other = parser.add_argument_group('Other options')
    g_other.add_argument('--version', action='version', version=verstr)
    g_other.add_argument(
        '-v',
        '--verbose',
        dest='verbose_count',
        action='count',
        default=0,
        help='increases log verbosity for each occurrence, debug level is -vvv',
    )
    g_other.add_argument(
        '-w',
        '--work-dir',
        action='store',
        type=Path,
        default=Path('work').absolute(),
        help='path where intermediate results should be stored',
    )
    g_other.add_argument(
        '--clean-workdir',
        action='store_true',
        default=False,
        help=(
            'Clears working directory of contents. Use of this flag is not'
            'recommended when running concurrent processes of aslprep.'
        ),
    )
    g_other.add_argument(
        '--resource-monitor',
        action='store_true',
        default=False,
        help="enable Nipype's resource monitoring to keep track of memory and CPU usage",
    )
    g_other.add_argument(
        '--config-file',
        action='store',
        metavar='FILE',
        help=(
            'Use pre-generated configuration file. Values in file will be overridden '
            'by command-line arguments.'
        ),
    )
    g_other.add_argument(
        '--write-graph',
        action='store_true',
        default=False,
        help='Write workflow graph.',
    )
    g_other.add_argument(
        '--stop-on-first-crash',
        action='store_true',
        default=False,
        help='Force stopping on first crash, even if a work directory was specified.',
    )
    g_other.add_argument(
        '--notrack',
        action='store_true',
        default=False,
        help=(
            'Opt-out of sending tracking information of this run to '
            'the aslprep developers. This information helps to '
            'improve aslprep and provides an indicator of real '
            'world usage crucial for obtaining funding.'
        ),
    )
    g_other.add_argument(
        '--debug',
        action='store',
        nargs='+',
        choices=config.DEBUG_MODES + ('all',),
        help="Debug mode(s) to enable. 'all' is alias for all available modes.",
    )

    latest = check_latest()
    if latest is not None and currentv < latest:
        print(
            f"""\
You are using aslprep-{currentv}, and a newer version of aslprep is available: {latest}.
Please check out our documentation about how and when to upgrade:
https://aslprep.readthedocs.io/en/latest/faq.html#upgrading""",
            file=sys.stderr,
        )

    _blist = is_flagged()
    if _blist[0]:
        _reason = _blist[1] or 'unknown'
        print(
            f"""\
WARNING: Version {config.environment.version} of aslprep (current) has been FLAGGED
(reason: {_reason}).
That means some severe flaw was found in it and we strongly
discourage its usage.""",
            file=sys.stderr,
        )

    return parser


def parse_args(args=None, namespace=None):
    """Parse args and run further checks on the command line."""
    import logging

    from niworkflows.utils.spaces import Reference, SpatialReferences

    parser = _build_parser()
    opts = parser.parse_args(args, namespace)
    if opts.config_file:
        skip = {} if opts.reports_only else {'execution': ('run_uuid',)}
        config.load(opts.config_file, skip=skip, init=False)
        config.loggers.cli.info(f'Loaded previous configuration file {opts.config_file}')

    config.execution.log_level = int(max(25 - 5 * opts.verbose_count, logging.DEBUG))
    config.from_dict(vars(opts), init=['nipype'])

    # Consistency checks
    force_set = set(config.workflow.force)
    ignore_set = set(config.workflow.ignore)
    if {'bbr', 'no-bbr'} <= force_set:
        msg = (
            'Cannot force and disable boundary-based registration at the same time. '
            'Remove `bbr` or `no-bbr` from the `--force` options.'
        )
        raise ValueError(msg)
    if 'fmap-jacobian' in force_set & ignore_set:
        msg = (
            'Cannot force and ignore fieldmap Jacobian correction. '
            'Remove `fmap-jacobian` from either the `--force` or the `--ignore` option.'
        )
        raise ValueError(msg)

    # Initialize --output-spaces if not defined
    if config.execution.output_spaces is None:
        config.execution.output_spaces = SpatialReferences(
            [Reference('MNI152NLin2009cAsym', {'res': 'native'})]
        )

    # Retrieve logging level
    build_log = config.loggers.cli

    # Load base plugin_settings from file if --use-plugin
    if opts.use_plugin is not None:
        import yaml

        with open(opts.use_plugin) as f:
            plugin_settings = yaml.safe_load(f)
        _plugin = plugin_settings.get('plugin')
        if _plugin:
            config.nipype.plugin = _plugin
            config.nipype.plugin_args = plugin_settings.get('plugin_args', {})
            config.nipype.nprocs = opts.nprocs or config.nipype.plugin_args.get(
                'n_procs', config.nipype.nprocs
            )

    # Resource management options
    # Note that we're making strong assumptions about valid plugin args
    # This may need to be revisited if people try to use batch plugins
    if 1 < config.nipype.nprocs < config.nipype.omp_nthreads:
        build_log.warning(
            f'Per-process threads (--omp-nthreads={config.nipype.omp_nthreads}) exceed '
            f'total threads (--nthreads/--n_cpus={config.nipype.nprocs})'
        )

    # Inform the user about the risk of using brain-extracted images
    if config.workflow.skull_strip_t1w == 'auto':
        build_log.warning(
            """\
Option ``--skull-strip-t1w`` was set to 'auto'. A heuristic will be \
applied to determine whether the input T1w image(s) have already been skull-stripped.
If that were the case, brain extraction and INU correction will be skipped for those T1w \
inputs. Please, BEWARE OF THE RISKS TO THE CONSISTENCY of results when using varying \
processing workflows across participants. To determine whether a participant has been run \
through the shortcut pipeline (meaning, brain extraction was skipped), please check the \
citation boilerplate. When reporting results with varying pipelines, please make sure you \
mention this particular variant of aslprep listing the participants for which it was \
applied."""
        )

    bids_dir = config.execution.bids_dir
    output_dir = config.execution.output_dir
    work_dir = config.execution.work_dir
    version = config.environment.version

    if config.execution.fs_subjects_dir is None:
        config.execution.fs_subjects_dir = output_dir / 'sourcedata' / 'freesurfer'

    if config.execution.aslprep_dir is None:
        config.execution.aslprep_dir = output_dir

    # Wipe out existing work_dir
    if opts.clean_workdir and work_dir.exists():
        from niworkflows.utils.misc import clean_directory

        build_log.info(f'Clearing previous aslprep working directory: {work_dir}')
        if not clean_directory(work_dir):
            build_log.warning(f'Could not clear all contents of working directory: {work_dir}')

    # Update the config with an empty dict to trigger initialization of all config
    # sections (we used `init=False` above).
    # This must be done after cleaning the work directory, or we could delete an
    # open SQLite database
    config.from_dict({})

    # Ensure input and output folders are not the same
    if output_dir == bids_dir:
        rec_path = bids_dir / 'derivatives' / f'aslprep-{version.split("+")[0]}'
        parser.error(
            'The selected output folder is the same as the input BIDS folder. '
            f'Please modify the output path (suggestion: {rec_path}).'
        )

    if bids_dir in work_dir.parents:
        parser.error(
            'The selected working directory is a subdirectory of the input BIDS folder. '
            'Please modify the output path.'
        )

    # Validate inputs
    if not opts.skip_bids_validation:
        from fmriprep.utils.bids import validate_input_dir

        build_log.info(
            'Making sure the input data is BIDS compliant (warnings can be ignored in most cases).'
        )
        validate_input_dir(
            config.environment.exec_env,
            opts.bids_dir,
            opts.participant_label,
            need_T1w=not config.execution.derivatives,
        )

    # Setup directories
    config.execution.log_dir = config.execution.aslprep_dir / 'logs'
    # Check and create output and working directories
    config.execution.log_dir.mkdir(exist_ok=True, parents=True)
    output_dir.mkdir(exist_ok=True, parents=True)
    work_dir.mkdir(exist_ok=True, parents=True)

    # Force initialization of the BIDSLayout
    config.execution.init()
    all_subjects = config.execution.layout.get_subjects()
    if config.execution.participant_label is None:
        config.execution.participant_label = all_subjects

    participant_label = set(config.execution.participant_label)
    missing_subjects = participant_label - set(all_subjects)
    if missing_subjects:
        parser.error(
            'One or more participant labels were not found in the BIDS directory: '
            f'{", ".join(missing_subjects)}.'
        )

    config.execution.participant_label = sorted(participant_label)
    config.workflow.skull_strip_template = config.workflow.skull_strip_template[0]
