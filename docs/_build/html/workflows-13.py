from aslprep.workflows.bold import init_bold_surf_wf
wf = init_bold_surf_wf(
    mem_gb=1,
    surface_spaces=['fsnative', 'fsaverage5'],
    medial_surface_nan=False)