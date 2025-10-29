"""
The workflow builder factory method.

All the checks and the construction of the workflow are done
inside this function that has pickleable inputs and output
dictionary (``retval``) to allow isolation using a
``multiprocessing.Process`` that allows aslprep to enforce
a hard-limited memory-scope.

"""


def build_workflow(config_file, retval):
    """Create the Nipype Workflow that supports the whole execution graph."""
    from fmriprep.utils.bids import check_pipeline_version
    from fmriprep.utils.misc import fmt_subjects_sessions
    from niworkflows.utils.misc import check_valid_fs_license

    from aslprep import config, data
    from aslprep.reports.core import generate_reports
    from aslprep.utils.misc import check_deps
    from aslprep.workflows.base import init_aslprep_wf

    config.load(config_file)
    build_log = config.loggers.workflow

    version = config.environment.version

    retval['return_code'] = 1
    retval['workflow'] = None

    banner = [f'Running ASLPrep version {version}']
    notice_path = data.load.readable('NOTICE')
    if notice_path.exists():
        banner[0] += '\n'
        banner += [f'License NOTICE {"#" * 50}']
        banner += [f'ASLPrep {version}']
        banner += notice_path.read_text().splitlines(keepends=False)[1:]
        banner += ['#' * len(banner[1])]
    build_log.log(25, f'\n{" " * 9}'.join(banner))

    # warn if older results exist: check for dataset_description.json in output folder
    msg = check_pipeline_version(
        'ASLPrep',
        version,
        config.execution.aslprep_dir / 'dataset_description.json',
    )
    if msg is not None:
        build_log.warning(msg)

    # Called with reports only
    if config.execution.reports_only:
        from aslprep.data import load as load_data

        build_log.log(
            25,
            'Running --reports-only on %s',
            fmt_subjects_sessions(config.execution.processing_groups),
        )
        session_list = config.execution.session_label
        if not session_list:
            session_list = (
                config.execution.bids_filters.get('asl', {}).get('session')
                if config.execution.bids_filters
                else None
            )

        failed_reports = generate_reports(
            config.execution.participant_label,
            config.execution.aslprep_dir,
            config.execution.run_uuid,
            session_list=session_list,
            bootstrap_file=load_data('reports-spec.yml'),
        )
        if failed_reports:
            config.loggers.cli.error(
                'Report generation was not successful for the following participants : '
                f'{", ".join(failed_reports)}.'
            )

        retval['return_code'] = len(failed_reports)
        return retval

    # Build main workflow
    init_msg = [
        "Building ASLPrep's workflow:",
        f'BIDS dataset path: {config.execution.bids_dir}.',
        f'Participants and sessions: {fmt_subjects_sessions(config.execution.processing_groups)}.',
        f'Run identifier: {config.execution.run_uuid}.',
        f'Output spaces: {config.execution.output_spaces}.',
    ]

    if config.execution.derivatives:
        init_msg += [f'Searching for derivatives: {config.execution.derivatives}.']

    if config.execution.fs_subjects_dir:
        init_msg += [f"Pre-run FreeSurfer's SUBJECTS_DIR: {config.execution.fs_subjects_dir}."]

    build_log.log(25, f'\n{" " * 11}* '.join(init_msg))

    retval['workflow'] = init_aslprep_wf()

    # Check for FS license after building the workflow
    if not check_valid_fs_license():
        build_log.critical(
            """\
ERROR: a valid license file is required for FreeSurfer to run. ASLPrep looked for an existing \
license file at several paths, in this order: 1) command line argument ``--fs-license-file``; \
2) ``$FS_LICENSE`` environment variable; and 3) the ``$FREESURFER_HOME/license.txt`` path. Get it \
(for free) by registering at https://surfer.nmr.mgh.harvard.edu/registration.html"""
        )
        retval['return_code'] = 126  # 126 == Command invoked cannot execute.
        return retval

    # Check workflow for missing commands
    missing = check_deps(retval['workflow'])
    if missing:
        build_log.critical(
            'Cannot run ASLPrep. Missing dependencies:%s',
            '\n\t* '.join([''] + [f'{cmd} (Interface: {iface})' for iface, cmd in missing]),
        )
        retval['return_code'] = 127  # 127 == command not found.
        return retval

    config.to_filename(config_file)
    build_log.info(
        'ASLPrep workflow graph with %d nodes built successfully.',
        len(retval['workflow']._get_all_nodes()),
    )
    retval['return_code'] = 0
    return retval


def build_boilerplate(config_file, workflow):
    """Write boilerplate in an isolated process."""
    from aslprep import config

    config.load(config_file)
    logs_path = config.execution.aslprep_dir / 'logs'
    boilerplate = workflow.visit_desc()
    citation_files = {ext: logs_path / f'CITATION.{ext}' for ext in ('bib', 'tex', 'md', 'html')}

    if boilerplate:
        # Rename "BOLD" to "ASL" in the boilerplate
        boilerplate = boilerplate.replace('BOLD', 'ASL')

        # To please git-annex users and also to guarantee consistency
        # among different renderings of the same file, first remove any
        # existing one
        for citation_file in citation_files.values():
            try:
                citation_file.unlink()
            except FileNotFoundError:
                pass

    citation_files['md'].write_text(boilerplate)

    if not config.execution.md_only_boilerplate and citation_files['md'].exists():
        from shutil import copyfile
        from subprocess import CalledProcessError, TimeoutExpired, check_call

        from aslprep.data import load as load_data

        # Generate HTML file resolving citations
        cmd = [
            'pandoc',
            '-s',
            '--bibliography',
            str(load_data('boilerplate.bib')),
            '--citeproc',
            '--metadata',
            'pagetitle="ASLPrep citation boilerplate"',
            str(citation_files['md']),
            '-o',
            str(citation_files['html']),
        ]

        config.loggers.cli.info('Generating an HTML version of the citation boilerplate...')
        try:
            check_call(cmd, timeout=10)
        except (FileNotFoundError, CalledProcessError, TimeoutExpired):
            config.loggers.cli.warning('Could not generate CITATION.html file:\n%s', ' '.join(cmd))

        # Generate LaTex file resolving citations
        cmd = [
            'pandoc',
            '-s',
            '--bibliography',
            str(load_data('boilerplate.bib')),
            '--natbib',
            str(citation_files['md']),
            '-o',
            str(citation_files['tex']),
        ]
        config.loggers.cli.info('Generating a LaTeX version of the citation boilerplate...')
        try:
            check_call(cmd, timeout=10)
        except (FileNotFoundError, CalledProcessError, TimeoutExpired):
            config.loggers.cli.warning('Could not generate CITATION.tex file:\n%s', ' '.join(cmd))
        else:
            copyfile(load_data('boilerplate.bib'), citation_files['bib'])
