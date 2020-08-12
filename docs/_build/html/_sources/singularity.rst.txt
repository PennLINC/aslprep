.. include:: links.rst


.. _run_singularity:

Running *ASLPrep* via Singularity containers
=============================================
Preparing a Singularity image
-----------------------------
**Singularity version >= 2.5**.
If the version of Singularity installed on your :abbr:`HPC (High-Performance Computing)`
system is modern enough you can create Singularity image directly on the system.
This is as simple as: ::

    $ singularity build /my_images/aslprep-<version>.simg docker://pennlinc/aslprep:<version>

where ``<version>`` should be replaced with the desired version of *ASLPrep* that you
want to download.

**Singularity version < 2.5**.
In this case, start with a machine (e.g., your personal computer) with Docker installed.
Use `docker2singularity <https://github.com/singularityware/docker2singularity>`_ to
create a singularity image.
You will need an active internet connection and some time. ::

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v D:\host\path\where\to\output\singularity\image:/output \
        singularityware/docker2singularity \
        pennlinc/aslprep:<version>

Where ``<version>`` should be replaced with the desired version of *ASLPrep* that you want
to download.

Beware of the back slashes, expected for Windows systems.
For \*nix users the command translates as follows: ::

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /absolute/path/to/output/folder:/output \
        singularityware/docker2singularity \
        pennlinc/aslprep:<version>


Transfer the resulting Singularity image to the HPC, for example, using ``scp``. ::

    $ scp pennlinc_aslprep*.img user@hcpserver.edu:/my_images

Running a Singularity Image
---------------------------
If the data to be preprocessed is also on the HPC, you are ready to run *ASLPrep*. ::

    $ singularity run --cleanenv aslprep.simg \
        path/to/data/dir path/to/output/dir \
        participant \
        --participant-label label

Handling environment variables
------------------------------
Singularity by default `exposes all environment variables from the host inside
the container <https://github.com/singularityware/singularity/issues/445>`__.
Because of this, your host libraries (e.g., nipype_ or a Python 2.7 environment)
could be accidentally used instead of the ones inside the container.
To avoid such a situation, we recommend using the ``--cleanenv`` argument in
all scenarios. For example: ::

    $ singularity run --cleanenv aslprep.simg \
      /work/04168/asdf/lonestar/ $WORK/lonestar/output \
      participant \
      --participant-label 387 --nthreads 16 -w $WORK/lonestar/work \
      --omp-nthreads 16


Alternatively, conflicts might be preempted and some problems mitigated by
unsetting potentially problematic settings, such as the ``PYTHONPATH`` variable,
before running: ::

    $ unset PYTHONPATH; singularity run aslprep.simg \
      /work/04168/asdf/lonestar/ $WORK/lonestar/output \
      participant \
      --participant-label 387 --nthreads 16 -w $WORK/lonestar/work \
      --omp-nthreads 16

It is possible to define environment variables scoped within the container by
using the ``SINGULARITYENV_*`` magic, in combination with ``--cleanenv``.
For example, we can set the FreeSurfer license variable 
as follows: ::

    $ export SINGULARITYENV_FS_LICENSE=$HOME/.freesurfer.txt
    $ singularity exec --cleanenv aslprep.simg env | grep FS_LICENSE
    FS_LICENSE=/home/users/oesteban/.freesurfer.txt

As we can see, the export in the first line tells Singularity to set a
corresponding environment variable of the same name after dropping the
prefix ``SINGULARITYENV_``.

Accessing the host's filesystem
-------------------------------
Depending on how Singularity is configured on your cluster it might or might not
automatically bind (mount or expose) host's folders to the container (e.g., ``/scratch``,
or ``$HOME``).
This is particularly relevant because, *if you can't run Singularity in privileged
mode* (which is almost certainly true in all the scenarios), **Singularity containers
are read only**.
This is to say that you won't be able to write *anything* unless Singularity can
access the host's filesystem in write mode.

By default, Singularity automatically binds (mounts) the user's *home* directory and
a *scratch* directory.
In addition, Singularity generally allows binding the necessary folders with
the ``-B <host_folder>:<container_folder>[:<permissions>]`` Singularity argument.
For example: ::

    $ singularity run --cleanenv -B /work:/work aslprep.simg \
      /work/my_dataset/ /work/my_dataset/derivatives/aslprep \
      participant \
      --participant-label 387 --nthreads 16 \
      --omp-nthreads 16

.. warning::

   If your Singularity installation doesn't allow you to bind non-existent bind points,
   you'll get an error saying ``WARNING: Skipping user bind, non existent bind point
   (directory) in container``.
   In this scenario, you can either try to bind things onto some other bind point you
   know it exists in the image or rebuild your singularity image with ``docker2singularity``
   as follows:
   ::

     $ docker run --privileged -ti --rm -v /var/run/docker.sock:/var/run/docker.sock \
              -v $PWD:/output singularityware/docker2singularity \
              -m "/gpfs /scratch /work /share /lscratch /opt/templateflow"

   In the example above, the following bind points are created: ``/gpfs``, ``/scratch``,
   ``/work``, ``/share``, ``/opt/templateflow``.

.. note::

   One great feature of containers is their confinement or isolation from the host
   system.
   Binding mount points breaks this principle, as the container has now access to
   create changes in the host.
   Therefore, it is generally recommended to use binding scarcely and granting
   very limited access to the minimum necessary resources.

**Relevant aspects of the** ``$HOME`` **directory within the container**.
By default, Singularity will bind the user's ``$HOME`` directory in the host
into the ``/home/$USER`` (or equivalent) in the container.
Most of the times, it will also redefine the ``$HOME`` environment variable and
update it to point to the corresponding mount point in ``/home/$USER``.
However, these defaults can be overwritten in your system.
It is recommended to check your settings with your system's administrators.
If your Singularity installation allows it, you can workaround the ``$HOME``
specification combining the bind mounts argument (``-B``) with the home overwrite
argument (``--home``) as follows: ::

    $ singularity run -B $HOME:/home/aslprep --home /home/aslprep \
          --cleanenv aslprep.simg <aslprep arguments>


.. _singularity_tf:

*TemplateFlow* and Singularity
------------------------------
:ref:`TemplateFlow` is a helper tool that allows *fMRIPrep* (or any other neuroimaging workflow)
to programmatically access a repository of standard neuroimaging templates.
In other words, *TemplateFlow* allows *fMRIPrep* to dynamically change the templates that
are used, e.g., in the atlas-based brain extraction step or spatial normalization.

Default settings in the Singularity image should get along with the Singularity
installation of your system.
However, deviations from the default configurations of your installation may break
this compatibility.
A particularly problematic case arises when the home directory is mounted in the
container, but the ``$HOME`` environment variable is not correspondingly updated.
Typically, you will experience errors like ``OSError: [Errno 30] Read-only file system``
or ``FileNotFoundError: [Errno 2] No such file or directory: '/home/aslprep/.cache'``.

If it is not explicitly forbidden in your installation, the first attempt to overcome this
issue is manually setting the ``$HOME`` directory as follows: ::

    $ singularity run --home $HOME --cleanenv aslprep.simg <aslprep arguments>

If the user's home directory is not automatically bound, then the second step would include
manually binding it as in the section above: ::

    $ singularity run -B $HOME:/home/aslprep --home /home/aslprep \
          --cleanenv aslprep.simg <aslprep arguments>

Finally, if the ``--home`` argument cannot be used, you'll need to provide the container with
writable filesystems where *TemplateFlow*'s files can be downloaded.
In addition, you will need to indicate *ASLPrep* to update the default paths with the new mount
points setting the ``SINGULARITYENV_TEMPLATEFLOW_HOME`` variable. ::

    $ export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow  # Tell ASLPrep the mount point
    $ singularity run -B <writable-path-on-host>:/opt/templateflow \
          --cleanenv aslprep.simg <aslprep arguments>

Internet access problems
------------------------
We have identified several conditions in which running *ASLPrep* might fail because
of spotty or impossible access to Internet.

If your compute node cannot have access to Internet, then you'll need to make sure
you run *ASLPrep* with the ``--notrack`` argument and pull down from TemplateFlow
all the resources that will be necessary.

With templateflow :: 
   $ export TEMPLATEFLOW_HOME=/path/to/keep/templateflow
   $ python -m pip install -U templateflow  # Install the client
   $ python
   >>> import templateflow.api
   >>> templateflow.api.TF_S3_ROOT = 'http://templateflow.s3.amazonaws.com'
   >>> api.get(‘MNI152NLin6Asym’)

Finally, run the singularity image binding the appropriate folder:
::

  $ export SINGULARITYENV_TEMPLATEFLOW_HOME=/templateflow
  $ singularity run -B ${TEMPLATEFLOW_HOME:-$HOME/.cache/templateflow}:/templateflow \
        --cleanenv aslprep.simg <aslprep  arguments>



Running Singularity on a SLURM system
-------------------------------------
An example of ``sbatch`` script to run *ASLPrep* on a SLURM system [#1]_ is given `below <singularity.html#sbatch-slurm>`__.
The submission script will generate one task per subject using a *job array*.
Submission is as easy as:
::

  $ export STUDY=/path/to/some/folder
  $ sbatch --array=1-$(( $( wc -l $STUDY/data/participants.tsv | cut -f1 -d' ' ) - 1 )) sbatch.slurm


.. literalinclude:: _static/sbatch.slurm
   :language: bash
   :name: sbatch.slurm
   :caption: **sbatch.slurm**:


.. [#1]  assuming that *job arrays* and Singularity are available
