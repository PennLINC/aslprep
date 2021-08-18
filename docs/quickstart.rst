.. include:: links.rst

---------------------
Quick Start Tutorial
---------------------

There are many to run *ASLPrep* for running but most have sensible defaults and don’t need to be changed. 
This page describes basic steps to run  *ASLPrep*. ASLPrep requires valid BIDS datasets and this is example of BIDS dataset in 
`openneuro <https://openneuro.org/datasets/ds000240/versions/2.0.0>`_.

Suppose the following data is available in the BIDS input:

sub-01
├── anat
│   ├── sub-01_T1w.json
│   └── sub-01_T1w.nii.gz
└── perf
    ├── sub-01_asl.json
    ├── sub-01_asl.nii.gz
    └── sub-01_aslcontext.tsv


One way to process these data would be to call *ASLprep*  like this:
aslprep  \
  /path/to/inputs /path/to/outputs participant  \
  --fs-license-file /path/to/license.txt


It can be run with  docker  :ref:`run_docker` or singularity image :ref:`run_singularity`.
