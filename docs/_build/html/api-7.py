from aslprep.workflows.bold.registration import init_bold_t1_trans_wf
wf = init_bold_t1_trans_wf(freesurfer=True,
                           mem_gb=3,
                           omp_nthreads=1)