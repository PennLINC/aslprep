.. include:: links.rst


.. _run_singularity:

Running *ASLPrep* via Singularity containers
=============================================
Preparing a Singularity Image
-----------------------------
**Singularity version >= 2.5**:
If the version of Singularity installed on your :abbr:`HPC (High-Performance Computing)`
system is modern enough you can create a Singularity image directly on the system
using the following command: ::

    $ singularity build aslprep-<version>.simg docker://pennlinc/aslprep:<version>

where ``<version>`` should be replaced with the desired version of *ASLPrep* that you
want to download.

Running a Singularity Image
---------------------------
If the data to be preprocessed is also on the HPC or a personal computer,
you are ready to run *ASLPrep*. ::

    $ singularity run --cleanenv aslprep.simg \
        path/to/data/dir path/to/output/dir \
        participant \
        --participant-label label

Handling environment variables
------------------------------
By default, Singularity interacts with all environment variables from the host.
The host libraries could accidentally conflict with singularity variables.
To avoid such a situation, it is recommended that you sue the ``--cleanenv or -e``
flag.
For instance: ::

    $ singularity run --cleanenv aslprep.simg \
      /work/789/asdf/ $WORK/output \
      participant \
      --participant-label 01

**Relevant aspects of the** ``$HOME`` **directory within the container**.
By default, Singularity will bind the user's ``$HOME`` directory on the host
into the ``/home/$USER`` directory (or equivalent) in the container.
Most of the times, it will also redefine the ``$HOME`` environment variable and
update it to point to the corresponding mount point in ``/home/$USER``.
However, these defaults can be overwritten in your system.
It is recommended that you check your settings with your system's administrators.
If your Singularity installation allows it, you can work around the ``$HOME``
specification combining the bind mounts argument (``-B``) with the home overwrite
argument (``--home``) as follows: ::

    $ singularity run -B $HOME:/home/aslprep --home /home/aslprep \
          --cleanenv aslprep.simg <aslprep arguments>
