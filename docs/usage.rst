.. include:: links.rst

###########
Usage Notes
###########


*****************************
Execution and the BIDS format
*****************************

The main input to *ASLPREP* is the path of the dataset that needs processing.

The input dataset is required to be in valid :abbr:`BIDS (Brain Imaging Data Structure)` format,
and it must include at least one T1w structural image.
We highly recommend that you validate your dataset with the free, online
`BIDS Validator <http://bids-standard.github.io/bids-validator/>`_.

.. important::
    Please note that ASL data in BIDS datasets should
    `already be scaled <https://bids-specification.readthedocs.io/en/v1.8.0/
    04-modality-specific-files/01-magnetic-resonance-imaging-data.html#scaling>`_.

    What this means is that the M0 scans in your dataset should preferably be scaled before running
    ASLPrep.

    Please see
    `the BIDS starter kit <https://bids-standard.github.io/bids-starter-kit/tutorials/asl.html>`_
    for information about converting ASL data to BIDS.

    If your data are not already scaled, you should use the ``--m0_scale`` parameter when running
    ASLPrep.

The exact command to run *ASLPrep* depends on the :doc:`installation` method.
The common parts of the command follow the `BIDS-Apps
<https://github.com/BIDS-Apps>`_ definition.
For example: ::

    aslprep data/bids_root/ out/ participant -w work/


**********************
Command-Line Arguments
**********************

.. argparse::
   :ref: aslprep.cli.parser._build_parser
   :prog: aslprep
   :nodefault:
   :nodefaultconst:


***************************************
Running *ASLPrep* via Docker containers
***************************************

For every new version of *ASLPrep* that is released, a corresponding Docker
image is generated.

In order to run *ASLPrep* Docker images, the Docker Engine must be installed.

If you have used *ASLPrep* via Docker in the past, you might need to pull down a
more recent version of the image::

    $ docker pull pennlinc/aslprep:<latest-version>

*ASLPrep* can be run interacting directly with the Docker Engine via the ``docker run``
command line, or through a lightweight wrapper that was created for convenience.


Running *ASLPrep* directly interacting with the Docker Engine
=============================================================

**Running containers as a user**.

In order to run docker smoothly, it is best to prevent permissions issues
associated with the root file system. Running docker as user on the host is to
ensure the ownership of files written during the container execution.

*ASLPrep* requires a significant amount of memory, typicaly around 12GB per subject.
If using docker desktop, you can set this in preferences. You can also set it on the command line.

A ``docker`` container can be created using the following command::

    $ docker run -ti -m 12GB --rm \
        -v path/to/data:/data:ro \
        -v path/to/output:/out \
        pennlinc/aslprep:<latest-version> \
        /data /out/out \
        participant

For example::

    $ docker run -ti -m 12GB --rm \
        -v $HOME/ds000240:/data:ro \
        -v $HOME/ds000240-results:/out:rw \
        -v $HOME/tmp/ds000240-workdir:/work \
        -v ${FREESURFER_HOME}:/fs \
        pennlinc/aslprep:<latest-version> \
        /data /out/aslprep-<latest-version> \
        participant \
        --participant-label '01'
        --fs-license-file ${FREESURFER_HOME}/license.txt
        -w /work

See :ref:`usage` for more information.


********************************************
Running *ASLPrep* via Singularity containers
********************************************


Preparing a Singularity Image
=============================

**Singularity version >= 2.5**:
If the version of Singularity installed on your :abbr:`HPC (High-Performance Computing)`
system is modern enough you can create a Singularity image directly on the system
using the following command: ::

    $ singularity build aslprep-<version>.simg docker://pennlinc/aslprep:<version>

where ``<version>`` should be replaced with the desired version of *ASLPrep* that you
want to download.


Running a Singularity Image
===========================

If the data to be preprocessed is also on the HPC or a personal computer,
you are ready to run *ASLPrep*. ::

    $ singularity run --cleanenv aslprep.simg \
        path/to/data/dir path/to/output/dir \
        participant \
        --participant-label label


Handling environment variables
==============================

By default, Singularity interacts with all environment variables from the host.
The host libraries could accidentally conflict with singularity variables.
To avoid such a situation, it is recommended that you sue the ``--cleanenv or -e``
flag.
For instance: ::

    $ singularity run --cleanenv aslprep.simg \
        /work/789/asdf/ $WORK/output \
        participant \
        --participant-label 01

**Relevant aspects of the** ``$HOME`` **directory within the container**.
By default, Singularity will bind the user's ``$HOME`` directory on the host
into the ``/home/$USER`` directory (or equivalent) in the container.
Most of the times, it will also redefine the ``$HOME`` environment variable and
update it to point to the corresponding mount point in ``/home/$USER``.
However, these defaults can be overwritten in your system.
It is recommended that you check your settings with your system's administrators.
If your Singularity installation allows it, you can work around the ``$HOME``
specification combining the bind mounts argument (``-B``) with the home overwrite
argument (``--home``) as follows: ::

    $ singularity run -B $HOME:/home/aslprep --home /home/aslprep \
        --cleanenv aslprep.simg <aslprep arguments>


**********************
The FreeSurfer license
**********************

*ASLPRep* uses FreeSurfer tools, which require a license to run.

To obtain a FreeSurfer license, simply register for free at
https://surfer.nmr.mgh.harvard.edu/registration.html.

When using manually-prepared environments or singularity, FreeSurfer will search
for a license key file first using the ``$FS_LICENSE`` environment variable and then
in the default path to the license key file (``$FREESURFER_HOME/license.txt``).
If using the ``--cleanenv`` flag and ``$FS_LICENSE`` is set, use ``--fs-license-file $FS_LICENSE``
to pass the license file location to *ASLPrep*.

It is possible to run the docker container pointing the image to a local path
where a valid license file is stored.
For example, if the license is stored in the ``$HOME/.licenses/freesurfer/license.txt``
file on the host system: ::

    $ docker run -ti --rm \
        -v $HOME/fullds005:/data:ro \
        -v $HOME/dockerout:/out \
        -v $HOME/.licenses/freesurfer/license.txt:/opt/freesurfer/license.txt \
        pennlinc/aslprep:latest \
        /data /out/out \
        participant \
        --ignore fieldmaps


***************
Troubleshooting
***************

Logs and crashfiles are written to the
``<output dir>/aslprep/sub-<participant_label>/log`` directory.
Information on how to customize and understand these files can be found on the
`nipype debugging <http://nipype.readthedocs.io/en/latest/users/debug.html>`_
page.


*************************
Support and communication
*************************

The documentation of this project is found here: https://aslprep.readthedocs.io.

All bugs, concerns and enhancement requests for this software can be submitted here:
https://github.com/PennLINC/aslprep/issues.

If you have a question about using ``aslprep``,
please create a new topic on `NeuroStars <https://neurostars.org>`_ with
`the "Software Support" category and the "aslprep" tag
<https://neurostars.org/tags/c/software-support/234/aslprep>`_.
The ``aslprep`` developers follow NeuroStars, and will be able to answer your question there.
