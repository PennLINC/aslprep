"""Tests for aslprep.reports."""

from unittest.mock import patch

from aslprep.reports.core import run_reports


def test_run_reports_error_handling(tmp_path):
    """Report generation errors are written to the log directory.

    Reproduces nipreps/fmriprep#3636.
    """
    with patch('aslprep.reports.core.Report') as MockReport:
        MockReport.return_value.generate_report.side_effect = Exception('Test Exception')

        res = run_reports(
            output_dir=str(tmp_path),
            subject_label='01',
            run_uuid='test_uuid',
            errorname='report.err',
        )

        assert res == '01'
        error_file = tmp_path / 'logs' / 'report.err'
        assert error_file.is_file()

        content = error_file.read_text()
        assert 'Traceback' in content
        assert 'Test Exception' in content


def test_generate_reports_per_subject(tmp_path, monkeypatch):
    """Report specs and session lists are chosen per subject (nipreps/fmriprep#3409)."""
    from unittest.mock import MagicMock

    from aslprep import config
    from aslprep.reports import core

    sessions = {'01': [], '02': ['a', 'b', 'c', 'd', 'e']}
    layout = MagicMock()
    layout.get_sessions.side_effect = lambda subject, **kwargs: sessions[subject]
    monkeypatch.setattr(config.execution, 'layout', layout)
    monkeypatch.setattr(config.execution, 'aggr_ses_reports', 4)
    monkeypatch.setattr(config.execution, 'bids_filters', None)

    calls = []
    monkeypatch.setattr(
        core,
        'run_reports',
        lambda *args, **kwargs: calls.append(kwargs) or None,
    )
    core.generate_reports(['sub-01', '02'], tmp_path, 'uuid')

    out_files = [call['out_filename'] for call in calls]
    # sub-01 gets one aggregated report; sub-02 gets an anatomical report and session reports
    assert out_files[0] == 'report.html'
    assert out_files[1] == 'sub-02_anat.html'
    assert out_files[2:] == [f'sub-02_ses-{ses}_asl.html' for ses in sessions['02']]
    assert calls[1]['bootstrap_file'] != calls[0]['bootstrap_file']
