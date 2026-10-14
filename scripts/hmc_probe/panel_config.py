"""How each real-panel dataset is run: adapted copies, BIDS filters and extra options.

Decided in Task 0 of ``docs/plans/2026-10-13-hmc-probe.md`` ("Runs that will not get through
ASLPrep as curated"). The zoo itself is never modified: an adapted dataset is a copy made by
``adapt_panel.py``, with a note of every change.
"""

#: Datasets that are copied and changed before running. Each entry lists the changes, applied in
#: order by ``adapt_panel.py``.
ADAPTATIONS = {
    'ds004146': [
        {
            'op': 'remove_files',
            'glob': 'sub-*/ses-*/perf/*_m0scan.*',
            'why': 'M0Type is Included and a separate M0 is also associated; ASLPrep raises an '
            'error. The in-series M0 is kept.',
        },
    ],
    'ds006148': [
        {
            'op': 'drop_unpaired',
            'why': '11 controls and 10 labels; the unpaired final control is dropped.',
        },
    ],
    'ds000254': [
        {
            'op': 'drop_unpaired',
            'why': '44 controls and 45 labels; the unpaired final label is dropped.',
        },
    ],
    'kashyap': [
        {
            'op': 'unlink_m0',
            'glob': 'sub-*/ses-*/perf/*_dir-AP_*m0scan.json',
            'why': 'Two separate M0 candidates per run; the one with the same phase-encoding '
            'direction as the ASL runs (dir-PA) is kept.',
        },
    ],
    'jsr-aslprep-debug': [
        {
            'op': 'unlink_m0',
            'glob': 'sub-*/ses-03/perf/*_acq-calShortTR_m0scan.json',
            'why': 'Two separate M0 candidates in session 3; the long-TR calibration scan is kept.',
        },
    ],
}

#: Datasets left out of the benchmark runs, with the reason. Their originals are still run as
#: malformed-input tests.
EXCLUDED = {
    'ds007908': 'Background suppression with M0Type Absent: ASLPrep refuses to calibrate from '
    'suppressed controls at workflow construction, so no outputs exist to score.',
}

#: Extra command-line arguments per dataset.
EXTRA_ARGS = {
    'ds004873': ['--skip-bids-validation'],
    'ds008020': ['--skull-strip-t1w', 'skip'],
}

#: BIDS filters per dataset (written to a ``--bids-filter-file``).
BIDS_FILTERS = {
    'ds008719': {'asl': {'echo': 1}},
    'shanmugan_hup6': {'t1w': {'reconstruction': None, 'run': 1}},
}

#: Datasets whose calibration is known to be wrong: no CBF-level metrics.
NO_CBF_METRICS = {
    'ds000240': 'M0 RepetitionTimePreparation is 0 (non-finite saturation correction, #654).',
    'ds005752': 'GE M0 scaling factor unknown.',
    'ds007522/sub-HML0202': 'M0 and ASL on different, undocumented intensity scales.',
    'ds004856/ses-wave2': 'Suspected control/label order error.',
    'ds004856/ses-wave3': 'Suspected control/label order error.',
}
