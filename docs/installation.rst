.. include:: links.rst

------------
Installation
------------

There are two ways to get *ASLPrep* installed:

* within a  `Manually Prepared Environment (Python 3.7+)`_ 

* using container technologies (RECOMMENDED), such as :ref:`run_docker`
  or :ref:`run_singularity`.

Once the  environment  is  set-up (first option above),
the next step is executing ``aslprep`` in the command-line.
The ``aslprep`` command-line options are documented in the :ref:`usage`
section.::

  $ aslprep  <input_bids_path> <derivatives_path> <analysis_level> <named_options>

For a docker  or singularity image, then the command-line will include binding configurations
and options. The container execution followed by the ``aslprep`` command-line options
as if you were running it on a normal installation.
The command-line structure above is then modified as follows:
::

  $ <container_command_and_options> <container_image> \
       <input_bids_path> <derivatives_path> <analysis_level> <ASLPrep_named_options>


Manually Prepared Environment (Python 3.7+)
===========================================

.. warning::

   This method is not recommended! Please check out container alternatives
   in :ref:`run_docker`, and :ref:`run_singularity`.

Make sure all of *ASLPRep*'s `External Dependencies`_ are installed.
These tools must be installed and their binaries available in the
system's ``$PATH``.
A relatively interpretable description of how your environment can be set up
is found in the `Dockerfile <https://github.com/pennlinc/aslprep/blob/master/Dockerfile>`_.
As an additional installation setting, FreeSurfer requires a license file (see :ref:`fs_license`).

In a functional Python 3.7 (or above) environment with ``pip`` installed,
*ASLPRep* can be installed using the habitual command ::

    $ python -m pip install aslprep

Check your installation with the ``--version`` argument ::

    $ aslprep --version


External Dependencies
---------------------

*ASLPRep* is written using Python 3.7 (or above), and is based on
nipype_.

*ASLPRep* requires some other neuroimaging software tools that are
not handled by Python's packaging system (Pypi) used to deploy
the ``ASLPrep`` package:

- FSL_ (version 5.0.9)
- ANTs_ (version 2.2.0 - NeuroDocker build)
- AFNI_ (version Debian-16.2.07)
- `C3D <https://sourceforge.net/projects/c3d/>`_ (version 1.0.0)
- FreeSurfer_ (version 6.0.1)
- `bids-validator <https://github.com/bids-standard/bids-validator>`_ (version 1.1.0)
- `connectome-workbench <https://www.humanconnectome.org/software/connectome-workbench>`_ (version Debian-1.3.2)

