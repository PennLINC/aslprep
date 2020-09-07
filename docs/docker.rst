.. include:: links.rst

.. _run_docker:

Running *ASLPrep* via Docker containers
========================================
For every new version of *ASLPrep* that is released, a corresponding Docker
image is generated.

In order to run *ASLPrep* Docker images, the Docker Engine must be installed
(see :ref:`installation_docker`).

If you have used *ASLPrep* via Docker in the past, you might need to pull down a
more recent version of the image: ::

    $ docker pull pennlinc/aslprep:<latest-version>

The *ASLPrep* can be run interacting directly with the Docker Engine via the `docker run`
command line, or through a lightweight wrapper that was created for convenience.

ref:`usage`),
automatically translating directories into Docker mount points.


Running *ASLPrep* directly interacting with the Docker Engine
--------------------------------------------------------------
**Running containers as a user**.

The best option to run  docker smoothly is to prevent permissions issues 
associated with root file system.Running the docker as user on the host  is to
ensure the ownership of files written during the container execution. it i

The ``docker`` can be invoke directly::

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
