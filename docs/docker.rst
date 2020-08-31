.. include:: links.rst

.. _run_docker:

Running *ASLPrep* via Docker containers
========================================
For every new version of *ASLPrep* that is released, a corresponding Docker
image is generated.
The Docker image *becomes* a container when the execution engine loads the
image and adds an extra layer that makes it *runnable*.
In order to run *ASLPrep* Docker images, the Docker Engine must be installed
(see :ref:`installation_docker`).

If you have used *ASLPrep* via Docker in the past, you might need to pull down a
more recent version of the image: ::

    $ docker pull pennlinc/aslprep:<latest-version>

You can run *ASLPrep* interacting directly with the Docker Engine via the `docker run`
command line, or you can use a lightweight wrapper we created for convenience.

ref:`usage`),
automatically translating directories into Docker mount points.


Running *ASLPrep* directly interacting with the Docker Engine
--------------------------------------------------------------
If you need a finer control over the container execution, or you feel comfortable
with the Docker Engine, you can avoid the extra software layer of the wrapper.

**Accessing filesystems in the host within the container**.
Containers are confined in a sandbox, so they can't access the host in any way, unless 
you explicitly prescribe acceptable accesses to the host.
The Docker Engine provides mounting filesystems into the container with the ``-v`` argument
and the following syntax: ``-v some/path/in/host:/absolute/path/within/container:ro``, where the trailing ``:ro`` specifies that the mount is read-only. The mount permissions modifiers can be omitted, which means the mount will have read-write
permissions.

In general, you'll want to provide at least two mount-points: one read-only
for the input data and one read-write to store the outputs.
Potentially, you'll want to provide one or two more mount-points: one for the working
directory, in case you need to debug some issue or re-use pre-cached result, and
a :ref:`TemplateFlow` folder to preempt the download of your favorite templates in every
run.

**Running containers as a user**.
By default, Docker will run the container as **root**.
Some share systems may limit this feature and only allow running containers as a user.
When the container is run as **root**, output written to filesystems mounted from the
host will have the user id ``1000`` by default.
In other words, you'll need to be able to run as **root** in the host to change permissions
or manage these files.

Alternatively, running as a user allows you to prevent these permissions issues.
It is possible to run as a user with the ``-u`` argument.
In general, we will want to use the same user ID as the running user in the host to
ensure the ownership of files written during the container execution.
Therefore, you will generally run the container with ``-u $( id -u )``.

You may also invoke ``docker`` directly::

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

Once the Docker Engine arguments are written, the remainder of the command line follows
the :ref:`usage`.
In other words, the first section of the command line is equivalent to the ``aslprep``
executable in a bare-metal installation: ::

    $ docker run -ti --rm \                      | These lines
        -v $HOME/ds000240:/data:ro \                | are equivalent
        -v $HOME/ds000240/derivatives:/out \        | to the aslprep
        -v $HOME/tmp/ds000240-workdir:/work \       | executable.
        pennlinc/aslprep:<latest-version> \  |
        \
        /data /out/aslprep-<latest-version> \   | These lines
        participant \                            | correspond to
        -w /work                                 | aslprep arguments.
