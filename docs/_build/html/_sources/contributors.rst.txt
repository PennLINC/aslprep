.. include:: links.rst

------------------------
Contributing to ASLPREP
------------------------

This document explains how to prepare a new development environment and
update an existing environment, as necessary.

Development in Docker is encouraged, for the sake of consistency and
portability.
By default, work should be built off of `pennlinc/aslprep:unstable
<https://hub.docker.com/r/pennlinc/aslprep/>`_, which tracks the ``master`` branch,
or ``pennlinc/aslprep:latest``, which tracks the latest release version (see the
installation_ guide for the basic procedure for running).


Adding dependencies
===================
New dependencies to be inserted into the Docker image will either be
Python or non-Python dependencies.
Python dependencies may be added in three places, depending on whether
the package is large or non-release versions are required.
The image `must be rebuilt <#rebuilding-docker-image>`_ after any
dependency changes.

Python dependencies should generally be included in the ``REQUIRES``
list in `aslprep/__about__.py
<https://github.com/pennlinc/aslprep/aslprep/__about__.py#L87-L107>`_.
If the latest version in `PyPI <https://pypi.org/>`_ is sufficient,
then no further action is required.


Code-Server Development Environment (Experimental)
==================================================
To get the best of working with containers and having an interactive
development environment, we have an experimental setup with `code-server
<https://github.com/cdr/code-server>`_.

.. Note::
    We have `a video walking through the process
    <https://youtu.be/bkZ-NyUaTvg>`_ if you want a visual guide.

1. Build the Docker image
~~~~~~~~~~~~~~~~~~~~~~~~~
We will use the ``Dockerfile_devel`` file to build 
our development docker image::

    $ cd $HOME/projects/aslprep
    $ docker build -t aslprep_devel -f Dockerfile_devel .

2. Run the Docker image
~~~~~~~~~~~~~~~~~~~~~~~
We can start a docker container using the image we built (``aslprep_devel``)::

    $ docker run -it -p 127.0.0.1:8445:8080 -v ${PWD}:/src/aslprep aslprep_devel:latest

.. Note::
    If you are using windows shell, ${PWD} may not be defined, instead use the absolute
    path to your aslprep directory.

3. Copy aslprep.egg-info into your aslprep directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``aslprep.egg-info`` makes the aslprep  package exacutable inside the docker container.
Open a terminal in vscode and type the following::

    $ cp -R /src/aslprep.egg-info /src/aslprep/


Code-Server Development Environment Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The editor is `vscode <https://code.visualstudio.com/docs>`_

- There are several preconfigured debugging tests under
  the debugging icon in the activity bar

  - see `vscode debugging python <https://code.visualstudio.com/docs/python/debugging>`_
    for details.

- The ``gitlens`` and ``python`` extensions are preinstalled to improve
  the development experience in vscode.


Adding new features to the citation boilerplate
===============================================

The citation boilerplate is built by adding two dunder attributes
of workflow objects: ``__desc__`` and ``__postdesc__``.
Once the full *aslprep* workflow is built, starting from the
outer workflow and visiting all sub-workflows in topological
order, all defined ``__desc__`` are appended to the citation
boilerplate before descending into sub-workflows.
Once all the sub-workflows of a given workflow have
been visited, then the ``__postdesc__`` attribute is appended
and the execution pops out to higher level workflows.
The dunder attributes are written in Markdown language, and may contain
references.
To add a reference, just add a new Bibtex entry to the references
database (``/aslprep/data/boilerplate.bib``).
You can then use the Bibtex handle within the Markdown text.
For example, if the Bibtex handle is ``myreference``, a citation
will be generated in Markdown language with ``@myreference``.
To generate citations with parenthesis and/or additional content,
brackets should be used: e.g., ``[see @myreference]`` will produce
a citation like *(see Doe J. et al 2018)*.


An example of how this works is shown here: ::

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
    Head-motion parameters with respect to the BOLD reference
    (transformation matrices, and six corresponding rotation and translation
    parameters) are estimated before any spatiotemporal filtering using
    `mcflirt` [FSL {fsl_ver}, @mcflirt].
    """.format(fsl_ver=fsl.Info().version() or '<ver>')
