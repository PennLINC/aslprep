.. include:: links.rst

========================
FAQ, Tips, and Tricks
========================

Should I run quality control on my images before running *ASLPrep*?
--------------------------------------------------------------------
Yes, you should do so before any processing/analysis takes place.

Oftentimes (more often than we would like), images have fatal artifacts and problems.

Some exclusion criteria for data quality should be pre-specified before QC and any screening
of the original data.
Those exclusion criteria must be designed in agreement with the goals and challenges of the
experimental design.
For instance, when running cortical thickness analyses, images should be excluded when they present even the most subtle ghosts or other artifacts that may introduce biases in surface
reconstruction.
However, if the same data was planned to be used as a reference for spatial
normalization, those artifacts should be noted, but may not call for exclusion of the data.

When using publicly available datasets, an additional concern is that images may have already gone through
some kind of preprocessing (see next question).


What if I find some images have undergone some pre-processing already (e.g., my T1w image is already skull-stripped)?
---------------------------------------------------------------------------------------------------------------------
These images imply an unknown level of preprocessing (e.g., was it already bias-field corrected?),
which makes it difficult to decide on best-practices for further processing.
Hence, supporting such images was considered very low priority for *ASLPrep*.


My *ASLPrep* run is hanging...
-------------------------------
When running on Linux platforms (or containerized environments, because they are built around
Ubuntu), there is a Python bug that affects *ASLPrep* and drives the Linux kernel to kill
processes as a response to running out of memory.
Depending on the process killed by the kernel, *ASLPrep* may crash with a ``BrokenProcessPool``
error, or hang indefinitely, depending on settings.
While we are working on finding a solution that does not run up against this bug, this may take some
time.
This can be most easily resolved by allocating more memory to the process, if possible.

Additionally, consider using the ``--low-mem`` flag, which will make some memory optimizations at the cost of disk space in the working directory.

ERROR: it appears that ``recon-all`` is already running
-------------------------------------------------------
When running FreeSurfer_'s ``recon-all``, an error may say *it appears it is already running*.
FreeSurfer creates files (called ``IsRunning.{rh,lh,lh+rh}``, under the ``scripts/`` folder)
to determine whether it is already executing ``recon-all`` on that particular subject
in another process, compute node, etc.
If a FreeSurfer execution terminates abruptly, those files are not wiped out, and therefore,
the next time you try to execute ``recon-all``, FreeSurfer *thinks* it is still running.
The output you get from *ASLPrep* will contain something like the following: ::

  RuntimeError: Command:
  recon-all -autorecon2-volonly -openmp 8 -subjid sub-020 -sd /outputs/freesurfer -nogcareg -nocanorm -nocareg -nonormalization2 -nomaskbfs -nosegmentation -nofill
  Standard output:
  Subject Stamp: freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.1-f53a55a
  Current Stamp: freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.1-f53a55a
  INFO: SUBJECTS_DIR is /outputs/freesurfer
  Actual FREESURFER_HOME /opt/freesurfer
  -rw-rw-r-- 1 11239 users 207798 Apr  1 16:19 /outputs/freesurfer/sub-020/scripts/recon-all.log
  Linux 62324c0da859 4.4.0-142-generic #168-Ubuntu SMP Wed Jan 16 21:00:45 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux

  ERROR: it appears that recon-all is already running
  for sub-020 based on the presence of /outputs/freesurfer/sub-020/scripts/IsRunning.lh+rh. It could
  also be that recon-all was running at one point but
  died in an unexpected way. If it is the case that there
  is a process running, you can kill it and start over or
  just let it run. If the process has died, you should type:

  rm /outputs/freesurfer/sub-020/scripts/IsRunning.lh+rh

  and re-run. Or you can add -no-isrunning to the recon-all
  command-line. The contents of this file are:
  ----------------------------------------------------------
  ------------------------------
  SUBJECT sub-020
  HEMI    lh rh
  DATE Fri Mar 22 20:33:09 UTC 2019
  USER root
  HOST 622795a21a5f
  PROCESSID 55530
  PROCESSOR x86_64
  OS Linux
  Linux 622795a21a5f 4.4.0-142-generic #168-Ubuntu SMP Wed Jan 16 21:00:45 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
  $Id: recon-all,v 1.580.2.16 2017/01/18 14:11:24 oesteban Exp $
  ----------------------------------------------------------
  Standard error:

  Return code: 1


As suggested by the ``recon-all`` output message, deleting these files will enable
FreeSurfer to execute ``recon-all`` again.
In general, please be cautious of deleting files and mindful as to why a file exists.



How much CPU time and RAM should I allocate for a typical *ASLPrep* run?
-----------------------------------------------------------------------
We recommend  running *ASLPrep* by processing  one subject per container instance. A typical preprocessing run
without surface processing with FreeSurfer can be completed in less than 1 hour with 4 CPUs or in about 30 hour with 16 CPUs.
More than 16 CPUs do not translate into faster processing times for a single subject. About 8GB of memory should be
available for a single subject preprocessing run.

.. _upgrading:

A new version of *ASLPrep* has been published, when should I upgrade?
----------------------------------------------------------------------
We follow a philosophy of releasing updates very often, although the pace is slowing down
with the maturation of the software.
It is very likely that your version gets outdated over the extent of your study.
If that is the case (an ongoing study), then we discourage changing versions.
In other words, **the whole dataset should be processed with the same version (and, if applicable, the 
same container build) of *ASLPrep*.**

On the other hand, if the project is about to start, then we strongly recommend
using the latest version of the tool.

In any case, if you find your release listed as *flagged* in `this file
of our repo <https://github.com/pennlinc/ASLPrep/blob/master/.versions.json>`__,
then please update as soon as possible.

I'm running *ASLPrep* via Singularity containers - how can I troubleshoot problems?
------------------------------------------------------------------------------------
We have extended `this documentation <singularity.html>`__ to cover some of the most
frequent issues other Singularity users have faced.
Generally, users have found it difficult to `get TemplateFlow and Singularity to work
together <singularity.html#singularity-tf>`__.

What is *TemplateFlow* for?
---------------------------
*TemplateFlow* enables *ASLPrep* to generate preprocessed outputs spatially normalized to
a number of different neuroimaging templates (e.g., MNI).
For further details, please check `its documentation section <spaces.html#templateflow>`__.

.. _tf_no_internet:

How do you use TemplateFlow in the absence of access to the Internet?
---------------------------------------------------------------------
This is a fairly common situation in :abbr:`HPC (high-performance computing)`
systems, where the so-called login nodes have access to the Internet but
compute nodes are isolated, or in PC/laptop enviroments if you are travelling.
*TemplateFlow* will require Internet access the first time it receives a
query for a template resource that has not been previously accessed.
If you know which templates you are planning to use, you could
pre-fetch them using the Python client.
To do so, follow the next steps.

  1. By default, a mirror of *TemplateFlow* to store the resources will be
     created in ``$HOME/.cache/templateflow``.
     You can modify such a configuration with the ``TEMPLATEFLOW_HOME``
     environment variable, e.g.::

       $ export TEMPLATEFLOW_HOME=$HOME/.templateflow

  2. Install the client within your favorite Python 3 environment (this can
     be done in your login-node, or in a host with Internet access,
     without need for Docker/Singularity)::

       $ python -m pip install -U templateflow

  3. Use the ``get()`` utility of the client to pull down all the templates you will
     want to use. For example::

       $ python -c "from templateflow.api import get; get(['MNI152NLin2009cAsym', 'MNI152NLin6Asym', 'OASIS30ANTs', 'MNIPediatricAsym', 'MNIInfant'])"

After getting the resources you need, you will just need to make sure your
runtime environment is able to access the filesystem, at the location of your
*TemplateFlow home* directory.
If you are a Singularity user, please check out :ref:`singularity_tf`.
