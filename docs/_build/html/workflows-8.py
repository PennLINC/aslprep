from aslprep.workflows.bold import init_bold_stc_wf
wf = init_bold_stc_wf(
    metadata={'RepetitionTime': 2.0,
              'SliceTiming': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]},
    )