#!/usr/bin/env python
"""
The aslprep on Docker wrapper


This is a lightweight Python wrapper to run aslprep.
Docker must be installed and running. This can be checked
running ::

  docker info


"""
import sys
import os
import re
import subprocess

__version__ = '99.99.99'
__copyright__ = 'Copyright 2020,'
__credits__ = ['wait']
__bugreports__ = 'https://github.com/pennlinc/issues'


MISSING = """
Image '{}' is missing
Would you like to download? [Y/n] """
PKG_PATH = '/usr/local/miniconda/lib/python3.7/site-packages'
TF_TEMPLATES = (
    'MNI152Lin',
    'MNI152NLin2009cAsym',
    'MNI152NLin6Asym',
    'MNI152NLin6Sym',
    'MNIInfant',
    'MNIPediatricAsym',
    'NKI',
    'OASIS30ANTs',
    'PNC',
    'UNCInfant',
    'fsLR',
    'fsaverage',
    'fsaverage5',
    'fsaverage6',
)
NONSTANDARD_REFERENCES = (
    'anat',
    'T1w',
    'run',
    'func',
    'sbref',
    'fsnative'
)

# Monkey-patch Py2 subprocess
if not hasattr(subprocess, 'DEVNULL'):
    subprocess.DEVNULL = -3

if not hasattr(subprocess, 'run'):
    # Reimplement minimal functionality for usage in this file
    def _run(args, stdout=None, stderr=None):
        from collections import namedtuple
        result = namedtuple('CompletedProcess', 'stdout stderr returncode')

        devnull = None
        if subprocess.DEVNULL in (stdout, stderr):
            devnull = open(os.devnull, 'r+')
            if stdout == subprocess.DEVNULL:
                stdout = devnull
            if stderr == subprocess.DEVNULL:
                stderr = devnull

        proc = subprocess.Popen(args, stdout=stdout, stderr=stderr)
        stdout, stderr = proc.communicate()
        res = result(stdout, stderr, proc.returncode)

        if devnull is not None:
            devnull.close()

        return res
    subprocess.run = _run


# De-fang Python 2's input - we don't eval user input
try:
    input = raw_input
except NameError:
    pass


def check_docker():
    """Verify that docker is installed and the user has permission to
    run docker images.

    Returns
    -------
    -1  Docker can't be found
     0  Docker found, but user can't connect to daemon
     1  Test run OK
     """
    try:
        ret = subprocess.run(['docker', 'version'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    except OSError as e:
        from errno import ENOENT
        if e.errno == ENOENT:
            return -1
        raise e
    if ret.stderr.startswith(b"Cannot connect to the Docker daemon."):
        return 0
    return 1


def check_image(image):
    """Check whether image is present on local system"""
    ret = subprocess.run(['docker', 'images', '-q', image],
                         stdout=subprocess.PIPE)
    return bool(ret.stdout)


def check_memory(image):
    """Check total memory from within a docker container"""
    ret = subprocess.run(['docker', 'run', '--rm', '--entrypoint=free',
                          image, '-m'],
                         stdout=subprocess.PIPE)
    if ret.returncode:
        return -1

    mem = [line.decode().split()[1]
           for line in ret.stdout.splitlines()
           if line.startswith(b'Mem:')][0]
    return int(mem)


def merge_help(wrapper_help, target_help):
    # Matches all flags with up to one nested square bracket
    opt_re = re.compile(r'(\[--?[\w-]+(?:[^\[\]]+(?:\[[^\[\]]+\])?)?\])')
    # Matches flag name only
    flag_re = re.compile(r'\[--?([\w-]+)[ \]]')

    # Normalize to Unix-style line breaks
    w_help = wrapper_help.rstrip().replace('\r', '')
    t_help = target_help.rstrip().replace('\r', '')

    w_usage, w_details = w_help.split('\n\n', 1)
    w_groups = w_details.split('\n\n')
    t_usage, t_details = t_help.split('\n\n', 1)
    t_groups = t_details.split('\n\n')

    w_posargs = w_usage.split('\n')[-1].lstrip()
    t_posargs = t_usage.split('\n')[-1].lstrip()

    w_options = opt_re.findall(w_usage)
    w_flags = sum(map(flag_re.findall, w_options), [])
    t_options = opt_re.findall(t_usage)
    t_flags = sum(map(flag_re.findall, t_options), [])

    # The following code makes this assumption
    assert w_flags[:2] == ['h', 'version']
    assert w_posargs.replace(']', '').replace('[', '') == t_posargs

    # Make sure we're not clobbering options we don't mean to
    overlap = set(w_flags).intersection(t_flags)
    expected_overlap = {
        'anat-derivatives',
        'fs-license-file',
        'fs-subjects-dir',
        'h',
        'use-plugin',
        'version',
        'w',
    }

    assert overlap == expected_overlap, "Clobbering options: {}".format(
        ', '.join(overlap - expected_overlap))

    sections = []

    # Construct usage
    start = w_usage[:w_usage.index(' [')]
    indent = ' ' * len(start)
    new_options = sum((
        w_options[:2],
        [opt for opt, flag in zip(t_options, t_flags) if flag not in overlap],
        w_options[2:]
    ), [])
    opt_line_length = 79 - len(start)
    length = 0
    opt_lines = [start]
    for opt in new_options:
        opt = ' ' + opt
        olen = len(opt)
        if length + olen <= opt_line_length:
            opt_lines[-1] += opt
            length += olen
        else:
            opt_lines.append(indent + opt)
            length = olen
    opt_lines.append(indent + ' ' + t_posargs)
    sections.append('\n'.join(opt_lines))

    # Use target description and positional args
    sections.extend(t_groups[:2])

    for line in t_groups[2].split('\n')[1:]:
        content = line.lstrip().split(',', 1)[0]
        if content[1:] not in overlap:
            w_groups[2] += '\n' + line

    sections.append(w_groups[2])

    # All remaining sections, show target then wrapper (skipping duplicates)
    sections.extend(t_groups[3:] + w_groups[6:])
    return '\n\n'.join(sections)


def is_in_directory(filepath, directory):
    return os.path.realpath(filepath).startswith(
        os.path.realpath(directory) + os.sep)


def get_parser():
    """Defines the command line interface of the wrapper"""
    import argparse

    class ToDict(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            d = {}
            for kv in values:
                k, v = kv.split("=")
                d[k] = os.path.abspath(v)
            setattr(namespace, self.dest, d)

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)

    # Standard ASLPREP arguments
    parser.add_argument('bids_dir', nargs='?', type=os.path.abspath,
                        default='')
    parser.add_argument('output_dir', nargs='?', type=os.path.abspath,
                        default='')
    parser.add_argument('analysis_level', nargs='?', choices=['participant'],
                        default='participant')

    parser.add_argument('-h', '--help', action='store_true',
                        help="show this help message and exit")
    parser.add_argument('--version', action='store_true',
                        help="show program's version number and exit")

    # Allow alternative images (semi-developer)
    parser.add_argument('-i', '--image', metavar='IMG', type=str,
                        default='pennlinc/aslprep{}'.format(__version__),
                        help='image name')

    # Options for mapping files and directories into container
    # Update `expected_overlap` variable in merge_help() when adding to this
    g_wrap = parser.add_argument_group(
        'Wrapper options',
        'Standard options that require mapping files into the container')
    g_wrap.add_argument('-w', '--work-dir', action='store', type=os.path.abspath,
                        help='path where intermediate results should be stored')
    g_wrap.add_argument(
        '--output-spaces', nargs="*",
        help="""\
Standard and non-standard spaces to resample anatomical and functional images to. \
Standard spaces may be specified by the form \
``<TEMPLATE>[:res-<resolution>][:cohort-<label>][...]``, where ``<TEMPLATE>`` is \
a keyword (valid keywords: %s) or path pointing to a user-supplied template, and \
may be followed by optional, colon-separated parameters. \
Non-standard spaces (valid keywords: %s) imply specific orientations and sampling \
grids. \
Important to note, the ``res-*`` modifier does not define the resolution used for \
the spatial normalization.""" % (', '.join('"%s"' % s for s in TF_TEMPLATES),
                                 ', '.join(NONSTANDARD_REFERENCES)))

    g_wrap.add_argument(
        '--fs-license-file', metavar='PATH', type=os.path.abspath,
        default=os.getenv('FS_LICENSE', None),
        help='Path to FreeSurfer license key file. Get it (for free) by registering'
             ' at https://surfer.nmr.mgh.harvard.edu/registration.html')
    g_wrap.add_argument(
        '--fs-subjects-dir', metavar='PATH', type=os.path.abspath,
        help='Path to existing FreeSurfer subjects directory to reuse. '
             '(default: OUTPUT_DIR/freesurfer)')
    g_wrap.add_argument(
        '--anat-derivatives', metavar='PATH', type=os.path.abspath,
        help='Path to existing sMRIPrep/aslprep-anatomical derivatives to fasttrack '
             'the anatomical workflow.')
    g_wrap.add_argument(
        '--use-plugin', metavar='PATH', action='store', default=None,
        type=os.path.abspath, help='nipype plugin configuration file')

    # Developer patch/shell options
    g_dev = parser.add_argument_group(
        'Developer options',
        'Tools for testing and debugging aslprep')
    g_dev.add_argument('--patch', nargs="+", metavar="PACKAGE=PATH", action=ToDict,
                       help='local repository to use within container')
    g_dev.add_argument('--shell', action='store_true',
                       help='open shell in image instead of running aslprep')
    g_dev.add_argument('--config', metavar='PATH', action='store',
                       type=os.path.abspath, help='Use custom nipype.cfg file')
    g_dev.add_argument('-e', '--env', action='append', nargs=2, metavar=('ENV_VAR', 'value'),
                       help='Set custom environment variable within container')
    g_dev.add_argument('-u', '--user', action='store',
                       help='Run container as a given user/uid. Additionally, group/gid can be'
                            'assigned, (i.e., --user <UID>:<GID>)')
    g_dev.add_argument('--network', action='store',
                       help='Run container with a different network driver '
                            '("none" to simulate no internet connection)')

    return parser


def main():
    """Entry point"""

    parser = get_parser()
    # Capture additional arguments to pass inside container
    opts, unknown_args = parser.parse_known_args()

    # Set help if no directories set
    if (opts.bids_dir, opts.output_dir, opts.version) == ('', '', False):
        opts.help = True

    # Stop if no docker / docker fails to run
    check = check_docker()
    if check < 1:
        if opts.version:
            print('aslprep wrapper {!s}'.format(__version__))
        if opts.help:
            parser.print_help()
        if check == -1:
            print("aslprep: Could not find docker command... Is it installed?")
        else:
            print("aslprep: Make sure you have permission to run 'docker'")
        return 1

    # For --help or --version, ask before downloading an image
    if not check_image(opts.image):
        resp = 'Y'
        if opts.version:
            print('aslprep wrapper {!s}'.format(__version__))
        if opts.help:
            parser.print_help()
        if opts.version or opts.help:
            try:
                resp = input(MISSING.format(opts.image))
            except KeyboardInterrupt:
                print()
                return 1
        if resp not in ('y', 'Y', ''):
            return 0
        print('Downloading. This may take a while...')

    # Warn on low memory allocation
    mem_total = check_memory(opts.image)
    if mem_total == -1:
        print('Could not detect memory capacity of Docker container.\n'
              'Do you have permission to run docker?')
        return 1
    if not (opts.help or opts.version or '--reports-only' in unknown_args) and mem_total < 8000:
        print('Warning: <8GB of RAM is available within your Docker '
              'environment.\nSome parts of aslprep may fail to complete.')
        if '--mem_mb' not in unknown_args:
            resp = 'N'
            try:
                resp = input('Continue anyway? [y/N]')
            except KeyboardInterrupt:
                print()
                return 1
            if resp not in ('y', 'Y', ''):
                return 0

    ret = subprocess.run(['docker', 'version', '--format', "{{.Server.Version}}"],
                         stdout=subprocess.PIPE)
    docker_version = ret.stdout.decode('ascii').strip()

    command = ['docker', 'run', '--rm', '-it', '-e',
               'DOCKER_VERSION_8395080871=%s' % docker_version]

    # Patch working repositories into installed package directories
    if opts.patch:
        for pkg, repo_path in opts.patch.items():
            command.extend(
                ['-v', '{}:{}/{}:ro'.format(repo_path, PKG_PATH, pkg)]
            )

    if opts.env:
        for envvar in opts.env:
            command.extend(['-e', '%s=%s' % tuple(envvar)])

    if opts.user:
        command.extend(['-u', opts.user])

    if opts.fs_license_file:
        command.extend([
            '-v', '{}:/opt/freesurfer/license.txt:ro'.format(
                opts.fs_license_file)])

    main_args = []
    if opts.bids_dir:
        command.extend(['-v', ':'.join((opts.bids_dir, '/data', 'ro'))])
        main_args.append('/data')
    if opts.output_dir:
        if not os.path.exists(opts.output_dir):
            # create it before docker does
            os.makedirs(opts.output_dir)
        command.extend(['-v', ':'.join((opts.output_dir, '/out'))])
        main_args.append('/out')
    main_args.append(opts.analysis_level)

    if opts.fs_subjects_dir:
        command.extend(['-v', '{}:/opt/subjects'.format(opts.fs_subjects_dir)])
        unknown_args.extend(['--fs-subjects-dir', '/opt/subjects'])

    if opts.anat_derivatives:
        command.extend(['-v', '{}:/opt/smriprep/subjects'.format(opts.anat_derivatives)])
        unknown_args.extend(['--anat-derivatives', '/opt/smriprep/subjects'])

    if opts.work_dir:
        command.extend(['-v', ':'.join((opts.work_dir, '/scratch'))])
        unknown_args.extend(['-w', '/scratch'])

    # Check that work_dir is not a child of bids_dir
    if opts.work_dir and opts.bids_dir:
        if is_in_directory(opts.work_dir, opts.bids_dir):
            print(
                'The selected working directory is a subdirectory of the input BIDS folder. '
                'Please modify the output path.')
            return 1

    if opts.config:
        command.extend(['-v', ':'.join((
            opts.config, '/home/aslprep/.nipype/nipype.cfg', 'ro'))])

    if opts.use_plugin:
        command.extend(['-v', ':'.join((opts.use_plugin, '/tmp/plugin.yml',
                                        'ro'))])
        unknown_args.extend(['--use-plugin', '/tmp/plugin.yml'])

    if opts.output_spaces:
        spaces = []
        for space in opts.output_spaces:
            if space.split(':')[0] not in (TF_TEMPLATES + NONSTANDARD_REFERENCES):
                tpl = os.path.basename(space)
                if not tpl.startswith('tpl-'):
                    raise RuntimeError("Custom template %s requires a `tpl-` prefix" % tpl)
                target = '/home/aslprep/.cache/templateflow/' + tpl
                command.extend(['-v', ':'.join((os.path.abspath(space), target, 'ro'))])
                spaces.append(tpl[4:])
            else:
                spaces.append(space)
        unknown_args.extend(['--output-spaces'] + spaces)

    if opts.shell:
        command.append('--entrypoint=bash')

    if opts.network:
        command.append('--network=' + opts.network)

    command.append(opts.image)

    # Override help and version to describe underlying program
    # Respects '-i' flag, so will retrieve information from any image
    if opts.help:
        command.append('-h')
        targethelp = subprocess.check_output(command).decode()
        print(merge_help(parser.format_help(), targethelp))
        return 0
    elif opts.version:
        # Get version to be run and exit
        command.append('--version')
        ret = subprocess.run(command)
        return ret.returncode

    if not opts.shell:
        command.extend(main_args)
        command.extend(unknown_args)

    print("RUNNING: " + ' '.join(command))
    ret = subprocess.run(command)
    if ret.returncode:
        print("aslprep: Please report errors to {}".format(__bugreports__))
    return ret.returncode


if __name__ == '__main__':
    sys.exit(main())
