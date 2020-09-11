.. include:: links.rst

.. _run_docker:

Running *ASLPrep* via Docker containers
========================================
For every new version of *ASLPrep* that is released, a corresponding Docker
image is generated.

In order to run *ASLPrep* Docker images, the Docker Engine must be installed.

If you have used *ASLPrep* via Docker in the past, you might need to pull down a
more recent version of the image: ::

    $ docker pull pennlinc/aslprep:<latest-version>

*ASLPrep* can be run interacting directly with the Docker Engine via the `docker run`
command line, or through a lightweight wrapper that was created for convenience.


Running *ASLPrep* directly interacting with the Docker Engine
--------------------------------------------------------------
**Running containers as a user**.

In order to run docker smoothly, it is best to prevent permissions issues
associated with the root file system. Running docker as user on the host is to
ensure the ownership of files written during the container execution.

A ``docker`` container can be created using the following command::

    $ docker run -ti --rm \
        -v path/to/data:/data:ro \
        -v path/to/output:/out \
        pennlinc/aslprep:<latest-version> \
        /data /out/out \
        participant

For example: ::

    $ docker run -ti --rm \
        -v $HOME/ds000240:/data:ro \
        -v $HOME/ds000240/derivatives:/out \
        -v $HOME/tmp/ds000240-workdir:/work \
        pennlinc/aslprep:<latest-version> \
        /data /out/aslprep-<latest-version> \
        participant \
        -w /work

See :ref:`usage` for more information.
