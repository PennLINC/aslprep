from aslprep.workflows.bold.confounds import init_ica_aroma_wf
wf = init_ica_aroma_wf(
    mem_gb=3,
    metadata={'RepetitionTime': 1.0},
    omp_nthreads=1)