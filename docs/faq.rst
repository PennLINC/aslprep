.. include:: links.rst

========================
FAQ, Tips, and Tricks
========================


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

