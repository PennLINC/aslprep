.. include:: links.rst

################################
Frequently Asked Questions (FAQ)
################################

.. contents::
   :local:
   :depth: 1


.. _supported_data:

*************************************
What kinds of ASL data are supported?
*************************************

*ASLPrep* reads ASL data organized according to the `BIDS`_ specification for ASL
(stored in the ``perf`` directory, with an ``aslcontext.tsv`` file describing each volume).
If your data are not in BIDS yet, `dcm2niix <https://github.com/rordenlab/dcm2niix>`_
writes most of the required metadata, and
`ExploreASL <https://github.com/ExploreASL/ExploreASL>`_ includes an import tool
that converts DICOMs to ASL-BIDS.

Supported data
==============

-  **Labeling schemes**: pseudo-continuous (PCASL), continuous (CASL),
   and pulsed (PASL) labeling.
   PASL data must use a bolus cut-off technique (QUIPSS, QUIPSS II, or Q2TIPS).
-  **Delays**: single-delay data (one post-labeling delay or inversion time),
   and multi-delay data (several post-labeling delays or inversion times),
   including multi-delay (P)CASL in which the labeling duration also varies.
   See :ref:`cbf_preproc` for the models.
-  **Readouts**: 2D multi-slice readouts (with slice timing, the delays are shifted slice by slice)
   and 3D readouts (e.g., 3D GRASE or 3D stack-of-spirals).
-  **Volume types**: control and label volumes, as well as scanner-generated
   difference (``deltam``) and CBF (``cbf``) volumes.
-  **M0 calibration**: a separate M0 scan, M0 volumes included in the ASL series,
   a single M0 estimate in the sidecar (``M0Estimate``),
   or, when there is no M0 at all, the control volumes
   (only if the data were acquired without background suppression).

Unsupported data
================

The following kinds of ASL data are *not* supported by *ASLPrep*.
Other tools that can process them are listed where we know of one.

.. warning::

   This list, and the tools listed for each kind of data, may be out of date.
   *ASLPrep* and the other tools add support for new kinds of ASL data over time,
   so check each tool's documentation before choosing one.

-  **Multi-echo (multi-TE) ASL**, which is used to estimate water exchange across
   the blood-brain barrier.
   *ASLPrep* does not model the echo time,
   so multi-echo data should not be processed with it.
   `ExploreASL <https://github.com/ExploreASL/ExploreASL>`_ supports multi-TE data
   (through FSL's FABBER), and the
   `FABBER ASL models <https://github.com/physimals/fabber_models_asl>`_
   include a multi-TE model that can be run directly.
-  **Time-encoded (Hadamard-encoded) PCASL** in its encoded form.
   *ASLPrep* does not decode the sub-boluses.
   `ExploreASL <https://github.com/ExploreASL/ExploreASL>`_ can decode and quantify time-encoded
   data, and `BASIL <https://asl-docs.readthedocs.io>`_ (``oxford_asl``) can decode it and
   fit the decoded data.
-  **Look-Locker readouts, such as QUASAR and Turbo-QUASAR**, which sample several inversion
   times after each label.
   BASIL provides dedicated tools for these:
   `QUASIL <https://asl-docs.readthedocs.io/en/latest/quasil.html>`_ for QUASAR and
   ``toast`` for Turbo-QUASAR.
-  **PASL without a bolus cut-off** (e.g., FAIR or EPISTAR without QUIPSS II or Q2TIPS),
   because the bolus duration is unknown.
   `BASIL <https://asl-docs.readthedocs.io>`_ can estimate the bolus duration from
   multi-delay PASL data.
-  **Velocity-selective ASL (VSASL)**.
   BIDS does not currently define a velocity-selective labeling type,
   and we are not aware of a general-purpose, open-source pipeline that quantifies VSASL data.
-  **Vessel-encoded (vessel-selective) ASL** for perfusion territory mapping.
   *ASLPrep* treats every labeling volume as non-selective.

If you have a kind of ASL data that is not listed here,
please `open an issue <https://github.com/PennLINC/aslprep/issues>`_ to ask about it.
If you would like *ASLPrep* to support a new kind of ASL data,
we will probably ask you to share example data,
since we need data of that kind to develop and test the new feature.


.. _poor_brain_mask:

***************************************
What can I do if my brain mask is poor?
***************************************

*ASLPrep* computes the ASL brain mask from the ASL reference image with SynthStrip by default,
which does not work well on all ASL data.
To use the method *ASLPrep* used before SynthStrip (FSL's BET combined with AFNI's
``3dAutomask``), run *ASLPrep* again with ``--skull-strip-asl legacy``.
This method has worked better on some datasets, such as HCP-Aging.
See :ref:`asl_skullstrip` for how each method works.

If your data use prescan normalization, ``--disable-n4`` may also help,
since N4 bias field correction tends to be unstable on such data.


.. _upgrading:

*********
Upgrading
*********

When *ASLPrep* starts, it checks whether a newer release is available and prints a message
like the following if one is::

   You are using aslprep-25.2.1, and a newer version of aslprep is available: 26.0.3.

The latest release is determined from the
`ASLPrep image tags on Docker Hub <https://hub.docker.com/r/pennlinc/aslprep/tags>`_,
and the result is cached in ``~/.cache/aslprep/latest`` for 14 days.
The check never prevents *ASLPrep* from running, and it is skipped silently if Docker Hub cannot be
reached.

Should I upgrade?
=================

New releases fix bugs and often change the outputs, so we recommend using the latest release
when you start processing a new dataset.
Check the `release notes <https://github.com/PennLINC/aslprep/releases>`_ to see what has
changed since the version you are using.

However, **do not mix versions within a study**.
Different versions of *ASLPrep* can produce different results from the same data,
so every participant in an analysis should be processed with the same version.
If you upgrade in the middle of a study, reprocess all of the participants with the new version.

How do I upgrade?
=================

*ASLPrep* is distributed as a container image, so upgrading means pulling the image for the
new version.
With Docker:

.. code-block:: bash

   docker pull pennlinc/aslprep:<version>

With Apptainer or Singularity:

.. code-block:: bash

   apptainer build aslprep-<version>.sif docker://pennlinc/aslprep:<version>

where ``<version>`` is the release you want (for example, ``26.0.3``).
We recommend pinning an explicit version rather than using the ``latest`` tag,
so that you know exactly which version produced your outputs.
The version is also recorded in ``dataset_description.json`` in the output directory.

When you run the new version, use new output and working directories.
Working directories in particular must not be reused across versions,
because cached results from the old version may be silently reused.

See :doc:`installation` for more details.

What does it mean if my version has been FLAGGED?
=================================================

If a serious problem is found in a release, that release is added to the list of flagged
versions in
`.versions.json <https://github.com/PennLINC/aslprep/blob/main/.versions.json>`_,
and *ASLPrep* prints a warning with the reason when you run it::

   WARNING: Version 26.0.0 of aslprep (current) has been FLAGGED
   (reason: Loading of derivative datasets broken.).

You should not use a flagged version.
Upgrade to the latest release and reprocess any data that the problem affects.
