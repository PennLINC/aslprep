#!/usr/bin/env python3
import logging
from pathlib import Path

logging.addLevelName(25, 'IMPORTANT')  # Add a new level between INFO and WARNING
logging.addLevelName(15, 'VERBOSE')  # Add a new level between INFO and DEBUG
logger = logging.getLogger('sdcflows')


def get_parser():
    """Define the command line interface"""
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from .. import __version__ as _vstr

    parser = ArgumentParser(description='SDC Workflows',
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        'bids_dir', action='store', type=Path,
        help='the root folder of a BIDS dataset')
    parser.add_argument('output_dir', action='store', type=Path,
                        help='the output path for the outcomes of preprocessing and visual '
                             'reports')
    parser.add_argument('analysis_level', choices=['participant', 'group'], nargs='+',
                        help='processing stage to be run, "participant" means individual analysis '
                             'and "group" is second level analysis.')
    # optional arguments
    parser.add_argument('--version', action='version', version='v{}'.format(_vstr))

    # Options that affect how pyBIDS is configured
    g_bids = parser.add_argument_group('Options for filtering BIDS queries')
    g_bids.add_argument('--participant-label', action='store', type=str,
                        nargs='*', dest='subject', help='process only particular subjects')
    g_bids.add_argument('--task', action='store', type=str, nargs='*',
                        help='select a specific task to be processed')
    g_bids.add_argument('--dir', action='store', type=str, nargs='*',
                        help='select a specific direction entity to be processed')
    g_bids.add_argument('--acq', action='store', type=str, nargs='*', dest='acquisition',
                        help='select a specific acquisition entity to be processed')
    g_bids.add_argument('--run', action='store', type=int, nargs='*',
                        help='select a specific run identifier to be processed')
    g_bids.add_argument('--suffix', action='store', type=str, nargs='*', default='bold',
                        help='select a specific run identifier to be processed')

    g_perfm = parser.add_argument_group('Options to handle performance')
    g_perfm.add_argument("-v", "--verbose", dest="verbose_count", action="count", default=0,
                         help="increases log verbosity for each occurence, debug level is -vvv")
    g_perfm.add_argument('--ncpus', '--nprocs', action='store', type=int,
                         help='maximum number of threads across all processes')
    g_perfm.add_argument('--nthreads', '--omp-nthreads', action='store', type=int,
                         help='maximum number of threads per-process')

    g_other = parser.add_argument_group('Other options')
    g_other.add_argument('-w', '--work-dir', action='store', type=Path,
                         help='path where intermediate results should be stored')

    return parser


def main():
    """Entry point"""
    from os import cpu_count
    from multiprocessing import set_start_method
    # from bids.layout import BIDSLayout
    from nipype import logging as nlogging
    set_start_method('forkserver')

    opts = get_parser().parse_args()

    # Retrieve logging level
    log_level = int(max(25 - 5 * opts.verbose_count, logging.DEBUG))
    # Set logging
    logger.setLevel(log_level)
    nlogging.getLogger('nipype.workflow').setLevel(log_level)
    nlogging.getLogger('nipype.interface').setLevel(log_level)
    nlogging.getLogger('nipype.utils').setLevel(log_level)

    # Resource management options
    plugin_settings = {
        'plugin': 'MultiProc',
        'plugin_args': {
            'n_procs': opts.ncpus,
            ''
            'raise_insufficient': False,
            'maxtasksperchild': 1,
        }
    }
    # Permit overriding plugin config with specific CLI options
    if not opts.ncpus or opts.ncpus < 1:
        plugin_settings['plugin_args']['n_procs'] = cpu_count()

    nthreads = opts.nthreads
    if not nthreads or nthreads < 1:
        nthreads = cpu_count()

    # output_dir = opts.output_dir.resolve()
    # bids_dir = opts.bids_dir or output_dir.parent

    # Get absolute path to BIDS directory
    # bids_dir = opts.bids_dir.resolve()
    # layout = BIDSLayout(str(bids_dir), validate=False, derivatives=str(output_dir))
    # query = {'suffix': opts.suffix, 'extension': ['.nii', '.nii.gz']}

    # for entity in ('subject', 'task', 'dir', 'acquisition', 'run'):
    #     arg = getattr(opts, entity, None)
    #     if arg is not None:
    #         query[entity] = arg


if __name__ == '__main__':
    raise RuntimeError("sdcflows/cli/run.py should not be run directly;\n"
                       "Please `pip install` sdcflows and use the `sdcflows` command")
