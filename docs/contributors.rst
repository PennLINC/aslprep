.. include:: links.rst

#########################
Contributing to *ASLPrep*
#########################

*ASLPrep* is a project of the *NiPreps* Community,
which provides detailed `contributing guidelines <https://www.nipreps.org/community/>`_.
Please review the *NiPreps* contributing guidelines before contributing to *ASLPrep*.

This document picks up where the *NiPreps* contributing guidelines leave off.
Specifically, it explains how to prepare a new development environment and update an existing
environment, as necessary.

Development in Docker is encouraged, for the sake of consistency and portability.
By default, work should be built off of
`pennlinc/aslprep:unstable <https://hub.docker.com/r/pennlinc/aslprep/>`_,
which tracks the ``main`` branch,
or ``pennlinc/aslprep:latest``,
which tracks the latest release
(see :doc:`installation` for the basic procedure for running).


***************************************************
Contributing to ASLPrep without adding dependencies
***************************************************

In the most common case, you will want to modify *ASLPrep*'s Python code without adding any
dependencies (Python or not) to the Docker image.
In this situation, you can use the ``unstable`` Docker image without having to build a new Docker
image yourself.

1. Pull the ``unstable`` Docker image.

   .. code-block:: bash

      docker pull pennlinc/aslprep:unstable

2. Fork the ASLPrep repository to your GitHub account.
   For more information on contributing via the fork-and-branch approach,
   see `GitHub's contributing guide
   <https://docs.github.com/en/get-started/quickstart/contributing-to-projects>`_.

3. Clone your forked repository to your local machine.

4. Create a branch to work on.
   **Make sure your branch is up to date with ASLPrep's ``main`` branch before making any
   changes!**

5. Make changes to the codebase that you want to try out.

6. Test our your changes by running the Docker container.
   The trick here is to mount your modified version of *ASLPrep* into the Docker container,
   overwriting the container's version.
   This way, your Docker container will run using your modified code,
   rather than the original version.

   You can do this by running the Docker image as described in :doc:`usage`,
   but adding in a mount point for your code:

   .. code-block:: bash

      docker run \
         -v /path/to/local/aslprep:/opt/conda/envs/aslprep/lib/python3.12/site-packages/aslprep \
         pennlinc/aslprep:unstable \
         ...  # see the usage documentation for info on what else to include in this command

7. Push your changes to GitHub.

8. Open a pull request to PennLINC/ASLPrep's ``main`` branch.
   Please follow `NiPreps contributing guidelines <https://www.nipreps.org/community/>`_
   when preparing a pull request.


Running tests
=============

While CircleCI will automatically run *ASLPrep*'s tests for any open PRs,
we strongly recommend running at least some tests locally, to make sure your proposed changes work.

*ASLPrep* has a file, ``aslprep/tests/run_local_tests.py``, that builds Docker ``run`` commands to
run selected tests.
Please use that script to run some tests locally before opening your PR.


********************************
Adding or modifying dependencies
********************************

If you think *ASLPrep* needs to use a library (Python or not) that is not installed in the Docker
image already, then you will need to build a new Docker image to test out your proposed changes.

*ASLPrep*'s Docker image is built from two Dockerfiles:

- ``Dockerfile.base`` contains non-Python runtime dependencies (FreeSurfer, AFNI, MSM, system
  libraries). It is rebuilt infrequently and published as a base image.
- ``Dockerfile`` uses `Pixi <https://pixi.sh>`_ to install the full Python environment
  (conda + PyPI dependencies) from ``pixi.lock`` and ``pyproject.toml``, then assembles the
  final image on top of the base.

Conda-level dependencies (FSL, ANTs, Connectome Workbench, TensorFlow, etc.) are declared in
``[tool.pixi.dependencies]`` in ``pyproject.toml``. Python dependencies are in
``[project.dependencies]``. Both are resolved together into ``pixi.lock``.

To add or modify a dependency, edit ``pyproject.toml`` and run ``pixi lock`` (on Linux) to
regenerate ``pixi.lock``, then rebuild the Docker image to verify::

  docker build --target aslprep -t pennlinc/aslprep:dev .

If you changed ``Dockerfile.base``, build the base image first (see ``.maint/INSTRUCTIONS.md``)::

  docker build -f Dockerfile.base -t pennlinc/aslprep_build:0.0.21 .

Once your change is working, open a pull request to the *ASLPrep* repo.
