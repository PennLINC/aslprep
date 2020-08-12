from aslprep.workflows.bold.registration import init_bbreg_wf
wf = init_bbreg_wf(use_bbr=True, bold2t1w_dof=9,
                   bold2t1w_init='register', omp_nthreads=1)