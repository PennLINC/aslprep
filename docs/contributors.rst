.. include:: links.rst

------------------------
Contributing to ASLPREP
------------------------

This document explains how to prepare a new development environment and
update an existing environment, as necessary.

Development in Docker is encouraged, for the sake of consistency and portability.
By default, work should be built off of `pennlinc/aslprep:unstable
<https://hub.docker.com/r/pennlinc/aslprep/>`_, which tracks the ``master`` branch,
or ``pennlinc/aslprep:latest``, which tracks the latest release version (see the
:ref:`installation` guide for the basic procedure for running).


Adding dependencies
===================
New dependencies to be inserted into the Docker image will either be
Python or non-Python dependencies.
Python dependencies may be added in three places, depending on whether
the package is large or if non-release versions are required.
The image `must be rebuilt <#rebuilding-docker-image>`_ after any
dependency changes.

Python dependencies should generally be included in the ``REQUIRES``
list in `aslprep/__about__.py
<https://github.com/pennlinc/aslprep/aslprep/__about__.py#L87-L107>`_.
If the latest version in `PyPI <https://pypi.org/>`_ is sufficient,
then no further action is required.


Code-Server Development Environment (Experimental)
==================================================

1. Build the Docker image
~~~~~~~~~~~~~~~~~~~~~~~~~
We will use the ``Dockerfile_devel`` file to build
our development docker image::

    $ mkdir $HOME/projects
    $ cd $HOME/projects/aslprep
    $ git clone https://github.com/PennLINC/aslprep.git
    $ docker build -t aslprep_devel -f Dockerfile_devel .

2. Run the Docker image
~~~~~~~~~~~~~~~~~~~~~~~
We can start a docker container using the image we built (``aslprep_devel``)::

    $ docker run -it -p 127.0.0.1:8445:8080 -v ${PWD}:/src/aslprep aslprep_devel:latest

.. Note::
    If you are using windows shell, ${PWD} may not be defined. Instead, use the absolute
    path to your *ASLPrep* directory.
