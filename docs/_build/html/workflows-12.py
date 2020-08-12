from aslprep.workflows.bold import init_bold_reg_wf
wf = init_bold_reg_wf(
    freesurfer=True,
    mem_gb=1,
    omp_nthreads=1,
    use_bbr=True,
    bold2t1w_dof=9,
    bold2t1w_init='register')