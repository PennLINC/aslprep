.. include:: links.rst

################################
Frequently Asked Questions (FAQ)
################################

.. contents::
   :local:
   :depth: 1


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
