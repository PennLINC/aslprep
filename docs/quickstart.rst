.. include:: links.rst

---------------------
Quick Start Tutorial
---------------------

This page describes the basic steps to run  *ASLPrep*. *ASLPrep* requires a valid BIDS dataset, like this example of BIDS dataset on 
`openneuro <https://openneuro.org/datasets/ds000240/versions/2.0.0>`_. The dataset can be downloaded from openneuro, either using 
command line or using GUI. With suscessful installation of `openneuro-cli <https://docs.openneuro.org/packages-openneuro-cli-readme>`_, 
the dataset can be downloaded using the command line::    

    openneuro download <accession number> <destination directory>

For example, to download a dataset with accession number ``ds000240``, type::
     
    openneuro download ds000240  /home/user/data
    
The BIDS dataset should include the following datatypes in order to run *ASLPrep*::

sub-01/anat/sub-01_T1w.nii.gz
sub-01/anat/sub-01_T1w.json
sub-01/perf/sub-01_asl.jnii.gz
sub-01/perf/sub-01_asl.json
sub-01/perf/sub-01_aslcontext.tsv


ASLPrep installation
---------------------
There are two ways to install *ASLPrep*:

1. Installation with ``pip``::

    python -m pip install aslprep
    
This method is not recommended, because it requires external dependencies to be installed. 


2. Installation through Docker/ Singularity
For every new version of *ASLPrep* that is released, a corresponding Docker
image is generated and pushed to `DockerHub <https://hub.docker.com/r/pennlinc/aslprep>`_. 
In order to run *ASLPrep* Docker images, the Docker Engine must be `installed <https://docs.docker.com/engine/install/>`_.

The docker image can be pulled from the ASLPrep DockerHub using the command line:: 

    docker pull pennlinc/aslprep:latest    

To use singularity,  a singularity image must be installed directly on the system
using the following command::

    singularity build aslprep.sif docker://pennlinc/aslprep:latest

This requires instalation of Singularity version >= 2.5

See  :ref:`run_docker`  and  :ref:`run_singularity` for more information. 

Running ASLPrep
----------------


If you are using a docker container, the following command can be called::

   docker run -ti --rm \                            # -ti: attach to the container interactively 
        -v path/to/data:/data:ro \                  # -v: mount the data directory to the container directory
        -v path/to/output:/out \                    # -v: mount the output directory to the container directory
        -v path/to/workdir:/work \                  # mount working directory
        -v ${FREESURFER_HOME}:/fs \                 # mount freesurfer home to get freesurfer license
        pennlinc/aslprep  \                         # the container name: aslprep
        /data /out/out \                            # the data and output directories
        participant                                 # analysis type: participant
        --participant-label 01                      # select participant 01
        --fs-license-file /path/to/license.txt      # setting freesurfer license file 
        -w /work                                    # setting working directory


For example, to run *ASLPrep* for one participant, run the following: ::

docker run -ti -m 12GB --rm -v $HOME/license.txt:/license/license.txt -v $HOME/ds000240:/data:ro -v $HOME/ds000240/derivatives:/out -v $HOME/tmp/ds000240-workdir:/work pennlinc/aslprep /data /out/aslprep participant --participant-label 01 --fs-license-file /license/license.txt -w /work

If you are using a singularity image, the following command can be called::

   singularity run --cleanenv aslprep.sif \       # run with singularity image  and clean environment  
        path/to/data/dir path/to/output/dir \     # the BIDS dataset directory and output directory 
        participant \                             # analysis  type by participant
        --participant-label label                 # select participant label eg, sub-01 or 01

Running *ASLPrep* will require a freesurfer license file, which can be requested `here <https://surfer.nmr.mgh.harvard.edu/registration.htm>`_.

For additional options, see usage notes >  :ref:`Usage`

ASLPrep outputs
---------------

 After a suscessful run, *ASLPrep* generates preprocessed ASL data, computed CBF maps, confound quality metrics, preprocessed 
 structural images, as well as one HTML report per subject that provides visual assessment of the 
 preprocessed data. See :ref:`Outputs` for more information. 

 








