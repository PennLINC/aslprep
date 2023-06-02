.. include:: links.rst

############
Installation
############

There are two ways to install *ASLPrep*:

* using `Container Technologies`_ (RECOMMENDED)
* within a `Manually Prepared Environment (Python 3.8+)`_, also known as *bare-metal installation*

*ASLPrep* is not just a Python package;
it also depends on a number of other neuroimaging libraries written in other languages.
As such, in order to install *ASLPrep* you also need to install these other dependencies
(with the right versions).
Container software makes this straightforward,
while the bare-metal approach is much more complicated and error-prone.

As such, we **strongly** recommend installing *ASLPrep* with Docker or Singularity.


.. _installation_container_technologies:

**********************
Container Technologies
**********************

*ASLPrep* is ideally run via a Docker or Singularity container.
If you are running *ASLPrep* locally, we suggest Docker_.
However, for security reasons, many :abbr:`HPCs (High-Performance Computing)` do not allow Docker
containers, but do allow Singularity_ containers.
The improved security for multi-tenant systems comes at the price of some limitations and extra
steps necessary for execution.


.. _installation_docker:

Docker Installation
===================

.. tip::

   *ASLPrep*'s image will take up ~30 GB of space on your machine.
   Please make sure to allocate Docker enough space to house the *ASLPrep* image before pulling it.

For every new version of *ASLPrep* that is released, a corresponding Docker image is generated and
deployed to DockerHub.

In order to run *ASLPrep* via Docker,
the Docker Engine `must be installed <https://docs.docker.com/engine/install/>`_.

If you have used *ASLPrep* via Docker in the past, you might need to pull down a more recent
version of the image:

.. code-block:: bash

   docker pull pennlinc/aslprep:<version>

The image can also be found here: https://registry.hub.docker.com/r/pennlinc/aslprep

*ASLPrep* can be run interacting directly with the Docker Engine via the ``docker run`` command,
or through a lightweight wrapper that was created for convenience.


.. _installation_singularity:

Singularity Installation
========================

**Singularity version >= 2.5**:
If the version of Singularity installed on your :abbr:`HPC (High-Performance Computing)` system is
modern enough, you can create a Singularity image directly on the system using the following
command:

.. code-block:: bash

   singularity build aslprep-<version>.simg docker://pennlinc/aslprep:<version>

where ``<version>`` should be replaced with the version of *ASLPrep* that you want to download.


.. _installation_manually_prepared_environment:

*******************************************
Manually Prepared Environment (Python 3.8+)
*******************************************

.. warning::

   This method is not recommended!
   Please use container alternatives described above instead.

ASLPrep requires some `External Dependencies`_.
These tools must be installed and their binaries available in the system's ``$PATH``.

On a functional Python 3.8 (or above) environment with ``pip`` installed,
ASLPrep can be installed using the command:

.. code-block:: bash

   pip install git+https://github.com/pennlinc/aslprep.git

Check your installation with the ``--version`` argument:

.. code-block:: bash

   aslprep --version


*********************
External Dependencies
*********************

*ASLPrep* is written using Python 3.8 (or above), is based on nipype_,
and requires some other neuroimaging software tools that are not handled by Python's packaging
system (PyPi) used to deploy the ASLPrep package:

-  FSL_ (version 6.0.2)
-  ANTs_ (version 2.2.0 - NeuroDocker build)
-  AFNI_ (version Debian-16.2.07)
-  `C3D <https://sourceforge.net/projects/c3d/>`_ (version 1.0.0)
-  FreeSurfer_ (version 6.0.1)
-  `bids-validator <https://github.com/bids-standard/bids-validator>`_ (version 1.1.0)
-  `connectome-workbench <https://www.humanconnectome.org/software/connectome-workbench>`_
   (version Debian-1.3.2)


*****************
Modifying ASLPrep
*****************

If you would like to make changes to ASLPrep's source code,
but do not wish to directly contribute to the package,
we recommend following the instructions in :doc:`contributors`.
