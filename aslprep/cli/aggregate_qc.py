# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Aggregate QC measures across all subjects in dataset."""

from pathlib import Path

from aslprep._warnings import logging


def _build_parser():
    """Build parser object."""
    from argparse import ArgumentParser, RawTextHelpFormatter
    from functools import partial

    def _path_exists(path, parser):
        """Ensure a given path exists."""
        if path is None or not Path(path).exists():
            raise parser.error(f'Path does not exist: <{path}>.')
        return Path(path).absolute()

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    PathExists = partial(_path_exists, parser=parser)

    parser.add_argument(
        'aslprep_dir',
        action='store',
        type=PathExists,
        help='Path to ASLPrep derivatives dataset.',
    )
    parser.add_argument(
        'output_dir',
        action='store',
        type=Path,
        help='Path to output directory where QC measures will be saved.',
    )
    parser.add_argument(
        'analysis_level',
        choices=['group'],
        help='processing stage to be run, only "group" in the case of this CLI.',
    )

    g_other = parser.add_argument_group('Other options')
    g_other.add_argument(
        '-w',
        '--work-dir',
        action='store',
        type=Path,
        default=Path('work').absolute(),
        help='path where intermediate results should be stored',
    )

    return parser


def parse_args(args=None, namespace=None):
    """Parse args and run further checks on the command line."""
    parser = _build_parser()
    opts = parser.parse_args(args, namespace)

    output_dir = opts.output_dir
    input_dir = opts.aslprep_dir
    work_dir = opts.work_dir

    output_dir.mkdir(exist_ok=True, parents=True)
    work_dir.mkdir(exist_ok=True, parents=True)

    log_dir = output_dir / 'logs'
    log_dir.mkdir(exist_ok=True, parents=True)

    build_log = logging.getLogger()

    # Wipe out existing work_dir
    if opts.clean_workdir and work_dir.exists():
        from niworkflows.utils.misc import clean_directory

        build_log.info(f'Clearing previous aslprep working directory: {work_dir}')
        if not clean_directory(work_dir):
            build_log.warning(f'Could not clear all contents of working directory: {work_dir}')

    return opts


def build_workflow(opts, retval):
    """Build the workflow."""
    from bids.layout import BIDSLayout
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from aslprep.data import load as load_data
    from aslprep.interfaces.bids import DerivativesDataSink
    from aslprep.interfaces.group import (
        AggregateCBFQC,
        MeanCBF,
        PlotAggregatedCBFQC,
        StatisticalMapRPT,
    )

    aslprep_config = load_data('aslprep_bids_config.json')
    layout = BIDSLayout(
        opts.aslprep_dir,
        validate=False,
        config=['bids', 'derivatives', aslprep_config],
    )
    qc_data = collect_aslprep_qc_derivatives(layout)

    workflow = Workflow(name='aslprep_group_wf')

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'cbf_qc',
                'cbf_img',
            ],
        ),
        name='inputnode',
    )
    inputnode.inputs.cbf_qc = qc_data['cbf_qc']
    inputnode.inputs.cbf_img = qc_data['cbf_img']

    # Aggregate QC metrics
    aggregate_cbf_qc = pe.Node(
        AggregateCBFQC(),
        name='aggregate_cbf_qc',
    )
    workflow.connect([(inputnode, aggregate_cbf_qc, [('cbf_qc', 'in_files')])])

    # Plot aggregated QC metrics
    plot_cbf_qc = pe.Node(
        PlotAggregatedCBFQC(),
        name='plot_cbf_qc',
    )
    workflow.connect([(aggregate_cbf_qc, plot_cbf_qc, [('out_file', 'in_file')])])

    # Save aggregated QC metrics
    ds_report_qei = pe.MapNode(
        DerivativesDataSink(
            base_directory=opts.output_dir,
            datatype='figures',
            suffix='cbf',
        ),
        iterfield=['in_file', 'desc'],
        name='ds_report_qei',
    )
    workflow.connect([
        (plot_cbf_qc, ds_report_qei, [
            ('out_files', 'in_file'),
            ('measures', 'desc'),
        ]),
    ])  # fmt:skip

    # Average CBF maps across subjects
    average_cbf_maps = pe.Node(
        MeanCBF(),
        name='average_cbf_maps',
    )
    workflow.connect([(inputnode, average_cbf_maps, [('cbf_img', 'in_files')])])

    # Plot average and SD CBF maps
    plot_mean_cbf = pe.Node(
        StatisticalMapRPT(cmap='viridis'),
        name='plot_mean_cbf',
    )
    workflow.connect([
        (inputnode, plot_mean_cbf, [('template', 'underlay')]),
        (average_cbf_maps, plot_mean_cbf, [('mean_map', 'overlay')]),
    ])  # fmt:skip

    ds_report_mean_cbf = pe.Node(
        DerivativesDataSink(
            base_directory=opts.output_dir,
            datatype='figures',
            desc='mean',
            suffix='cbf',
        ),
        name='ds_report_mean_cbf',
    )
    workflow.connect([(plot_mean_cbf, ds_report_mean_cbf, [('out_report', 'in_file')])])

    plot_sd_cbf = pe.Node(
        StatisticalMapRPT(cmap='viridis'),
        name='plot_sd_cbf',
    )
    workflow.connect([
        (inputnode, plot_sd_cbf, [('template', 'underlay')]),
        (average_cbf_maps, plot_sd_cbf, [('sd_map', 'overlay')]),
    ])  # fmt:skip

    ds_report_sd_cbf = pe.Node(
        DerivativesDataSink(
            base_directory=opts.output_dir,
            datatype='figures',
            desc='sd',
            suffix='cbf',
        ),
        name='ds_report_sd_cbf',
    )
    workflow.connect([(plot_sd_cbf, ds_report_sd_cbf, [('out_report', 'in_file')])])

    return workflow


def collect_aslprep_qc_derivatives(layout):
    """Collect ASLPrep QC derivatives."""
    qc_data = {}
    qc_files = layout.get(
        desc='qualitycontrol',
        suffix='cbf',
        extension='.tsv',
        return_type='file',
    )
    qc_data['cbf_qc'] = qc_files
    for qc_file in qc_files:
        # Collect associated files
        ...

    return qc_data


def build_boilerplate(opts, aslprep_wf):
    """Build the boilerplate."""
    pass


def main():
    """Entry point."""
    import gc
    import sys
    from multiprocessing import Manager, Process
    from os import EX_SOFTWARE
    from pathlib import Path

    from aslprep.utils.bids import write_bidsignore, write_derivative_description

    opts = parse_args()
    input_dir = opts.aslprep_dir
    output_dir = opts.output_dir
    work_dir = opts.work_dir

    logger = logging.getLogger()

    # CRITICAL Call build_workflow(config_file, retval) in a subprocess.
    # Because Python on Linux does not ever free virtual memory (VM), running the
    # workflow construction jailed within a process preempts excessive VM buildup.
    with Manager() as mgr:
        retval = mgr.dict()
        p = Process(target=build_workflow, args=(opts, retval))
        p.start()
        p.join()
        retval = dict(retval.items())  # Convert to base dictionary

        if p.exitcode:
            retval['return_code'] = p.exitcode

    exitcode = retval.get('return_code', 0)
    aslprep_wf = retval.get('workflow', None)

    exitcode = exitcode or (aslprep_wf is None) * EX_SOFTWARE
    if exitcode != 0:
        sys.exit(exitcode)

    # Generate boilerplate
    with Manager() as mgr:
        p = Process(target=build_boilerplate, args=(opts, aslprep_wf))
        p.start()
        p.join()

    # Clean up master process before running workflow, which may create forks
    gc.collect()

    logger.log(25, 'ASLPrep started!')
    errno = 1  # Default is error exit unless otherwise set
    try:
        aslprep_wf.run(**config.nipype.get_plugin())
    except Exception as e:
        logger.critical('ASLPrep-Group failed: %s', e)
        raise
    else:
        logger.log(25, 'ASLPrep-Group finished successfully!')

        # Bother users with the boilerplate only iff the workflow went okay.
        boiler_file = output_dir / 'logs' / 'CITATION.md'
        if boiler_file.exists():
            if config.environment.exec_env in (
                'singularity',
                'docker',
                'aslprep-docker',
            ):
                boiler_file = Path('<OUTPUT_PATH>') / boiler_file.relative_to(output_dir)
            logger.log(
                25,
                'Works derived from this ASLPrep execution should include the '
                f'boilerplate text found in {boiler_file}.',
            )
        errno = 0
    finally:
        from aslprep import data
        from aslprep.reports.core import generate_reports

        # Generate reports phase
        session_list = config.execution.get().get('bids_filters', {}).get('asl', {}).get('session')

        failed_reports = generate_reports(
            config.execution.participant_label,
            config.execution.aslprep_dir,
            config.execution.run_uuid,
            session_list=session_list,
            bootstrap_file=data.load('reports-spec.yml'),
        )
        write_derivative_description(input_dir, output_dir)
        write_bidsignore(output_dir)

        if failed_reports:
            msg = (
                'Report generation was not successful for the following participants '
                f': {", ".join(failed_reports)}.'
            )
            logger.error(msg)

        sys.exit(int((errno + len(failed_reports)) > 0))


if __name__ == '__main__':
    raise RuntimeError(
        'aslprep/cli/aggregate_qc.py should not be run directly;\n'
        'Please `pip install` aslprep and use the `aslprep-group` command'
    )
