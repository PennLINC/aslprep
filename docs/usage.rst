.. include:: links.rst

.. _Usage :

Usage Notes
===========

Execution and the BIDS format
-----------------------------
The main input to *ASLPREP* is the path of the dataset
that needs processing.

.. note::
     The official `BIDS Specification for ASL <https://docs.google.com/document/d/15tnn5F10KpgHypaQJNNGiNKsni9035GtDqJzWqkkP6c/edit>`_  is currently in its final stages of development.
     We have created a simple tool for converting ASL data into currently valid BIDS, available `at PennLINC <https://github.com/PennLINC/aslbids>`_.
     Note that some of the parameters necessary for running *ASLPrep* cannot be extracted directly from the DICOM header —
     to obtain these parameters, we recommend consulting the MRI physicist/technician at the scanner.
     The conversion of ASL into BIDS on `Flywheel <https://flywheel.io/>`_ using `fw-heudiconv <https://fw-heudiconv.readthedocs.io>`_
     is described `here at fw-heudiconv <https://fw-heudiconv.readthedocs.io/en/latest/tips.html#arterial-spin-labelling-data>`_.



The input dataset is required to be in valid :abbr:`BIDS (Brain Imaging Data
Structure)` format, and it must include at least one T1w structural image.
We highly recommend that you validate your dataset with the free, online
`BIDS Validator <http://bids-standard.github.io/bids-validator/>`_.

The exact command to run *ASLPrep* depends on the Installation_ method.
The common parts of the command follow the `BIDS-Apps
<https://github.com/BIDS-Apps>`_ definition.
For example: ::

    aslprep data/bids_root/ out/ participant -w work/



Command-Line Arguments
----------------------

.. argparse::
   :ref: aslprep.cli.parser._build_parser
   :prog: aslprep
   :nodefault:
   :nodefaultconst:

The FreeSurfer license
----------------------

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

Troubleshooting
---------------

Logs and crashfiles are written to the
``<output dir>/aslprep/sub-<participant_label>/log`` directory.
Information on how to customize and understand these files can be found on the
`nipype debugging <http://nipype.readthedocs.io/en/latest/users/debug.html>`_
page.

**Support and communication**.

The documentation of this project is found here: http://aslprep.readthedocs.org/en/latest/.

All bugs, concerns and enhancement requests for this software can be submitted here:
https://github.com/pennlinc/aslprep/issues.
