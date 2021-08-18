.. include:: links.rst

---------------------
Quick Start Tutorial
---------------------

There are many to run *ASLPrep* for running but most have sensible defaults and don’t need to be changed. 
This page describes basic steps to run  *ASLPrep*. ASLPrep requires valid BIDS datasets and this is example of BIDS dataset in 
`openneuro <https://openneuro.org/datasets/ds000240/versions/2.0.0>`_.

Suppose the following data is available in the BIDS input:: 

    sub-01/anat/sub-01_T1w.nii.gz
    sub-01/anat/sub-01_T1w.json
    sub-01/perf/sub-01_asl.jnii.gz
    sub-01/perf/sub-01_asl.json
    sub-01/perf/sub-01_aslcontext.tsv



ASLPrep installation
---------------------
There are two ways to install *ASLPrep*:

1. Installation through `pip`::
    
     pip install aslprep 
  
This method is recommended for users who are familiar with Python,have Python packages  and all dependencies installed. 
see :ref:`Installtation`.  Python 3.7 (or above) is required

2. Installation through Docker/ Singularity
For every new version of *ASLPrep* that is released, a corresponding Docker
image is generated. In order to run *ASLPrep* Docker images, the Docker Engine must be installed.

The docker image can be pull down as below:: 

    docker pull pennlinc/aslprep:latest 

For singularity image,  a singularity image directly on the system
using the following command::

    singularity build aslprep.sif docker://pennlinc/aslprep:latest

It requires instalation of Singularity version >= 2.5


See  :ref:`run_docker`  and  :ref:`run_singularity` for more information. 

Running ASLPrep
----------------

One way to process these data would be to call *ASLprep*  like this::

       aslprep  \
       /path/to/inputs /path/to/outputs participant  \
       --fs-license-file /path/to/license.txt \
       --participant_label sub-01  # for subject label sub-01

For any other options see the usage notes >  :ref:`Usage`

With a docker container, the follwoing command can be called::

   docker run -ti --rm \
        -v path/to/data:/data:ro \
        -v path/to/output:/out \
        pennlinc/aslprep  \
        /data /out/out \
        participant

For example: ::

    docker run -ti -m 12GB --rm \
        -v $HOME/ds000240:/data:ro \   # ds000240 is the BIDS dataset directory 
        -v $HOME/ds000240/derivatives:/out \
        -v $HOME/tmp/ds000240-workdir:/work \
        -v ${FREESURFER_HOME}:/fs \
        pennlinc/aslprep \
        /data /out/aslprep \
        participant \
        --participant-label 01
        --fs-license-file /fs/license.txt
        -w /work


With singularity image, the following command can be called::

   singularity run --cleanenv aslprep.sif \
        path/to/data/dir path/to/output/dir \
        participant \
        --participant-label label




ASLPrep outputs
---------------

 After suscessful run, *ASLPrep* generates outputs with  HTML report per subject that provide visual assessment of the 
 processed data. The outputs include preprocessed ASL data, computed CBF maps, cofound quality metrics, preprocessed 
 structural images, and other outputs. 

 








