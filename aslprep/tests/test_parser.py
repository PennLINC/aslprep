"""Test parser."""

from pathlib import Path

import pytest
from packaging.version import Version

from aslprep import config
from aslprep.cli import version as _version
from aslprep.cli.parser import _build_parser
from aslprep.tests.tests import reset_config

MIN_ARGS = ['data/', 'out/', 'participant']


@pytest.mark.parametrize(
    ('args', 'code'),
    [
        ([], 2),
        (MIN_ARGS, 2),  # bids_dir does not exist
        (MIN_ARGS + ['--fs-license-file'], 2),
        (MIN_ARGS + ['--fs-license-file', 'fslicense.txt'], 2),
    ],
)
def test_parser_errors(args, code):
    """Check behavior of the parser."""
    with pytest.raises(SystemExit) as error:
        _build_parser().parse_args(args)

    assert error.value.code == code


@pytest.mark.parametrize('args', [MIN_ARGS, MIN_ARGS + ['--fs-license-file']])
def test_parser_valid(tmp_path, args):
    """Check valid arguments."""
    datapath = tmp_path / 'data'
    datapath.mkdir(exist_ok=True)
    args[0] = str(datapath)

    if '--fs-license-file' in args:
        _fs_file = tmp_path / 'license.txt'
        _fs_file.write_text('')
        args.insert(args.index('--fs-license-file') + 1, str(_fs_file.absolute()))

    opts = _build_parser().parse_args(args)

    assert opts.bids_dir == datapath


@pytest.mark.parametrize(
    ('argval', 'gb'),
    [
        ('1G', 1),
        ('1GB', 1),
        ('1000', 1),  # Default units are MB
        ('32000', 32),  # Default units are MB
        ('4000', 4),  # Default units are MB
        ('1000M', 1),
        ('1000MB', 1),
        ('1T', 1000),
        ('1TB', 1000),
        (f'{1e6:.0f}K', 1),
        (f'{1e6:.0f}KB', 1),
        (f'{1e9:.0f}B', 1),
    ],
)
def test_memory_arg(tmp_path, argval, gb):
    """Check the correct parsing of the memory argument."""
    datapath = tmp_path / 'data'
    datapath.mkdir(exist_ok=True)
    _fs_file = tmp_path / 'license.txt'
    _fs_file.write_text('')

    args = MIN_ARGS + ['--fs-license-file', str(_fs_file)] + ['--mem', argval]
    opts = _build_parser().parse_args(args)

    assert opts.memory_gb == gb


@pytest.mark.parametrize(('current', 'latest'), [('1.0.0', '1.3.2'), ('1.3.2', '1.3.2')])
def test_get_parser_update(monkeypatch, capsys, current, latest):
    """Make sure the out-of-date banner is shown."""
    expectation = Version(current) < Version(latest)

    def _mock_check_latest(*args, **kwargs):
        return Version(latest)

    monkeypatch.setattr(config.environment, 'version', current)
    monkeypatch.setattr(_version, 'check_latest', _mock_check_latest)

    _build_parser()
    captured = capsys.readouterr().err

    msg = f"""\
You are using aslprep-{current}, and a newer version of aslprep is available: {latest}.
Please check out our documentation about how and when to upgrade:
https://aslprep.readthedocs.io/en/latest/faq.html#upgrading"""

    assert (msg in captured) is expectation


@pytest.mark.parametrize('flagged', [(True, None), (True, 'random reason'), (False, None)])
def test_get_parser_blacklist(monkeypatch, capsys, flagged):
    """Make sure the blacklisting banner is shown."""

    def _mock_is_bl(*args, **kwargs):
        return flagged

    monkeypatch.setattr(_version, 'is_flagged', _mock_is_bl)

    _build_parser()
    captured = capsys.readouterr().err

    assert ('FLAGGED' in captured) is flagged[0]
    if flagged[0]:
        assert (flagged[1] or 'reason: unknown') in captured


def test_reuse_config(tmp_path):
    """Check which settings are reused with ``--config-file``.

    Reproduces nipreps/fmriprep#3625.
    """
    from niworkflows.utils.testing import generate_bids_skeleton

    from aslprep.cli.parser import parse_args
    from aslprep.data import load as load_data

    reset_config()
    bids_dir = tmp_path / 'ds000240'
    generate_bids_skeleton(
        bids_dir,
        {'01': {'anat': {'suffix': 'T1w'}, 'perf': [{'suffix': 'asl'}]}},
    )
    # Avoid requiring validator installation
    cli_args = [
        str(bids_dir),
        str(tmp_path / 'out'),
        'participant',
        '--skip-bids-validation',
        '--skip-parcellation',
    ]

    parse_args(cli_args)
    default_config = config.get(flat=True)
    reset_config()

    # Simulate a configuration file written by a previous run with a different output directory
    config_text = Path(load_data('../tests/data/config.toml')).read_text()
    config_text = config_text.replace(
        '[execution]\n',
        f'[execution]\naslprep_dir = "{tmp_path / "old_out"}"\n',
    )
    config_file = tmp_path / 'config.toml'
    config_file.write_text(config_text)
    config_args = ['--config-file', str(config_file)]
    parse_args(cli_args + config_args)
    reused_config = config.get(flat=True)
    # Reusing the config will apply same values
    assert reused_config['execution.output_spaces'] != default_config['execution.output_spaces']
    assert reused_config['execution.output_spaces'] == (
        'asl T1w MNI152NLin2009cAsym:res-native fsaverage:den-10k fsaverage:den-30k'
    )
    # But some will still differ
    assert reused_config['execution.aslprep_dir'] == str(tmp_path / 'out')
    assert reused_config['execution.log_dir'] not in config_file.read_text()
    assert reused_config['execution.run_uuid'] not in config_file.read_text()
    reset_config()

    overridden_args = (
        cli_args + config_args + ['--output-spaces', 'MNI152NLin6Asym', '--force', 'bbr']
    )
    # set new output directory
    overridden_args[1] = str(tmp_path / 'out2')
    parse_args(overridden_args)
    overridden_config = config.get(flat=True)

    # Passed in argument will override
    assert overridden_config['execution.output_spaces'] == 'MNI152NLin6Asym:res-native'
    assert 'bbr' in overridden_config['workflow.force']

    # But some values will still differ
    for v in ('execution.run_uuid', 'execution.aslprep_dir'):
        assert reused_config[v] != overridden_config[v]
    reset_config()
