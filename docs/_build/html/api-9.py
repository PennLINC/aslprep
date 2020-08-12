from aslprep.workflows.bold.registration import init_fsl_bbr_wf
wf = init_fsl_bbr_wf(use_bbr=True, bold2t1w_dof=9, bold2t1w_init='register')