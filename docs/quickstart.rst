.. include:: links.rst

---------------------
Quick Start Tutorial
---------------------

This page describes the basic steps to run ``ASLPrep`` on a BIDS dataset. ``ASLPrep`` is containerized and available via ``pip``, and thus can be run in a variety of ways. Here, we provide the most basic and user friendly workflow.

*ASLPrep* requires a valid BIDS dataset, like this example of a BIDS dataset on `openneuro <https://openneuro.org/datasets/ds000240/versions/2.0.0/download>`_. Using Chrome, you can download the data via the browser. Note that you might have to create a new folder into which you can download the data. You can also acquire data using using ``datalad`` or ``aws``.

The BIDS dataset should include the following datatypes in order to run *ASLPrep*::

    sub-01/anat/sub-01_T1w.nii.gz
    sub-01/anat/sub-01_T1w.json
    sub-01/perf/sub-01_asl.jnii.gz
    sub-01/perf/sub-01_asl.json
    sub-01/perf/sub-01_aslcontext.tsv

ASLPrep installation
---------------------
There are two ways to install *ASLPrep*:

1. Installation through Docker or Singularity (**recommended**)

For every new version of ``ASLPrep`` that is released, a corresponding Docker image is generated and pushed to `DockerHub <https://hub.docker.com/r/pennlinc/aslprep>`_. In order to run *ASLPrep* Docker images, the Docker Engine must be `installed <https://docs.docker.com/engine/install/>`_.

We recommend using Docker or Singularity to run ASLPrep. The docker image can be pulled from the ASLPrep DockerHub using the command line::

    docker pull pennlinc/aslprep:latest

To use singularity, a singularity image must be installed directly on the system
using the following command::

    singularity build aslprep.sif docker://pennlinc/aslprep:latest

This requires installation of Singularity version >= 2.5

See  :ref:`run_docker`  and  :ref:`run_singularity` for more information.

2. Installation is available via ``pip``::

    python -m pip install aslprep

This method is not recommended, because it requires external dependencies to be installed.

Running ASLPrep
----------------

Running ``ASLPrep`` will require a freesurfer license file (you do not actually need Freesurfer, though), which can be requested `here <https://surfer.nmr.mgh.harvard.edu/registration.html>`_. Move this license to a folder in your ``$HOME`` directory (to find the path to your home directory in the terminal, echo $HOME) and call it **license**.

In the Docker desktop application, please select Preferences > Resources > Advanced and select at least 12GB for RAM. Restart Docker.

Move the data directory to your ``$HOME`` directory (again, to find this location out, run this in the terminal: echo $HOME). Make sure it is called **ds000240**.

The following command, which should run in about 8 hours, can be called for a single participant::

    docker run -ti -m 12GB --rm -v $HOME/license.txt:/license/license.txt -v $HOME/ds000240:/data:ro -v $HOME/ds000240-results:/out:rw -v $HOME/tmp/ds000240-workdir:/work:rw pennlinc/aslprep /data /out participant --participant-label 01 --fs-license-file /license/license.txt -w /work

Here is a breakdown of this command::

    docker run -ti -m 12GB --rm                         # attach to the container interactively
    -v $HOME/license.txt:/license/license.txt           # mount the freesurfer license directory
    -v $HOME/ds000240:/data:ro                          # mount the data directory to the container directory
    -v $HOME/ds000240-results:/out:rw                   # mount the output directory to the container directory
    -v $HOME/tmp/ds000240-workdir:/work                 # mount working directory
    pennlinc/aslprep                                    # the container name: aslprep
    /data                                               # the data directory
    /out/aslprep                                        # the output directory
    participant                                         # analysis type: participant
    --participant-label 01                              # select participant 01
    --fs-license-file /license/license.txt              # setting freesurfer license file
    -w /work                                            # setting working directory

For additional options, see usage notes >  :ref:`Usage`

ASLPrep outputs
---------------

 After a successful run, *ASLPrep* generates preprocessed ASL data, computed CBF maps, confound quality metrics, preprocessed
 structural images, as well as one HTML report per subject that provides visual assessment of the
 preprocessed data. See :ref:`Outputs` for more information.










