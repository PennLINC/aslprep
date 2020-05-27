from aslprep.workflows.bold import init_bold_hmc_wf
wf = init_bold_hmc_wf(
    mem_gb=1,
    omp_nthreads=1)