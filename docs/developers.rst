.. include:: links.rst

#######################
Developer documentation
#######################

**********************************
Testing against a simulated truth
**********************************

Besides unit tests and integration runs on real data, ASLPrep is tested on datasets simulated
by `aslscan <https://github.com/PennLINC/aslscan>`_, an ASL simulator that writes the ground
truth (perfusion, arrival time, partial-volume fractions, head motion) next to a BIDS dataset.
The design is in ``docs/specs/2026-10-08-aslscan-ground-truth-tests-design.md``.

Tiers
=====

Fast tier (runs with the unit tests)
   ASLPrep's ``ExtractCBF`` and ``ComputeCBF`` interfaces run directly on small simulated
   datasets (``fast_*`` recipes) and are compared with an independent implementation of the
   documented model (``aslprep/tests/truth_models.py``), to 0.1 % (median).
   Deliberately wrong expectations (label and control swapped, slice timing reversed, a
   shifted delay, a scaled M0, the labeling efficiency ignored) must fail the same bound, which
   shows the comparison is sensitive. ``test_aslscan_conformance.py``,
   ``test_aslscan_interpretation.py`` and ``test_aslscan_fixtures.py`` hold this tier.

Integration tier (one CircleCI job per recipe)
   ASLPrep runs end to end on a full-brain recipe, and ``aslprep/tests/truth_scoring.py``
   scores its outputs: CBF against the documented model (Tier A) and against the physical
   truth (Tier B), head-motion correction against the simulated trajectory, coregistration
   against a known anatomical offset, and the alignment of T1w- and MNI-space outputs.
   Each ``test_aslscan_<recipe>.py`` runs ASLPrep once and checks one metric per test item.

Why two tiers of CBF accuracy
=============================

aslscan simulates the general kinetic model, in which label that has reached tissue relaxes
with tissue T1. ASLPrep's single-delay equation (Alsop et al. 2015) assumes blood T1 throughout,
so even a perfect implementation recovers about 0.8x the true grey-matter CBF and 0.4x the
white-matter CBF on these data. Tier A therefore checks ASLPrep against what its documented
model should give; Tier B checks the physical accuracy relative to that expectation, never
against an absolute value.

Generating the fixtures
=======================

Fixture generation runs on Linux (or WSL) and needs ``git``, ``cargo`` (through ``rustup``)
and Python with ``numpy``, ``nibabel``, ``scipy`` and ``templateflow``.

.. code-block:: bash

   # build aslscan at the pinned revisions (writes a stamp next to the binary)
   python -m aslprep.tests.aslscan_fixtures --build-aslscan ~/aslscan-build

   # generate every recipe into a data directory, then verify it
   python -m aslprep.tests.aslscan_fixtures --generate aslprep/tests/test_data \
       --aslscan ~/aslscan-build/target/release/aslscan
   python -m aslprep.tests.aslscan_fixtures --verify aslprep/tests/test_data

   # voxel counts of the scoring masks
   python -m aslprep.tests.aslscan_fixtures --masks aslprep/tests/test_data

With ``ASLSCAN`` pointing at the built binary, tests generate a missing fixture on demand.
``ASLPREP_REQUIRE_FIXTURES=1`` (set in CI) turns a missing or stale fixture into an error.
A fixture is current when its manifest's digest matches its inputs: the pins, the builder
modules, its phantom and its own recipe, so editing one recipe leaves the others valid.

Running the tests locally
=========================

Build the test image from your checkout once (only ``Dockerfile`` or ``pixi.lock`` changes
need a rebuild; the checkout is mounted into the container):

.. code-block:: bash

   docker build --target test -t pennlinc/aslprep:test .

Then run a marker the way CircleCI does:

.. code-block:: bash

   ASLPREP_REQUIRE_FIXTURES=1 python aslprep/tests/run_local_tests.py \
       -m aslscan_pcasl1pld --data-dir aslprep/tests/test_data

Outputs, including ``truth_score.json``, go to ``aslprep/tests/pytests/out/<test name>/``.

Reading a failure
=================

Each scored item fails with the quantity, its value, the bound and the reason for the bound,
for example::

   native.tier_a.median_dev = 0.08 ratio, expected <= 0.05 (preprocessing (HMC, M0
   registration) leaves the median within 5 % of the model)

``truth_score.json`` holds every metric, including report-only ones (for example the
per-voxel spread that motion correction adds to delta-M, and the bias near the brain-mask
edge), and is stored as a CircleCI artifact.

Bounds and expected failures
============================

Ceilings live in ``aslprep/tests/truth_bounds.py``. They are fixed from what each method
should achieve and are never fitted to ASLPrep's output; a ceiling is loosened only in a
commit that explains why the method cannot meet it. A missing metric or ceiling fails.

Known disagreements are strict expected failures whose size is asserted first, so a change in
either direction fails: for example the single-delay Q2TIPS quantification, and the absence of
any background-suppression loss when the sidecar has ``LabelingEfficiency`` (#706).

When simulating background suppression, set aslscan's ``[background_suppression]
inversion_efficiency`` to 0.975. aslscan scales the label difference by (1 - 2 x efficiency)
per pulse, so 0.975 keeps 95 % of the ASL signal per pulse, the white-paper value
(Alsop et al. 2015) that ASLPrep assumes. aslscan's default of 0.95 keeps only 90 %.

Adding a recipe
===============

1. Create ``aslprep/tests/data/aslscan/<name>/`` with ``asl.json``, ``aslcontext.tsv``,
   ``overlay.toml`` (``seed``, ``[acquisition] matrix`` and ``noise_variance`` are required),
   ``recipe.toml`` and ``SOURCE.md``. Recipes are discovered from their directories.
2. Regenerate the spec file and commit it with the recipe::

      python -m aslprep.tests.aslscan_fixtures --spec > .circleci/aslscan_fixtures.txt

3. Generate the fixture, check its mask populations, and add a ``test_aslscan_<name>.py``
   module, a pytest marker in ``pyproject.toml``, and a CircleCI matrix entry.

Updating aslscan or a phantom
=============================

Change the pins (``ASLSCAN_REF``, ``MRSIM_ACQ_REF``, ``RUST_TOOLCHAIN``) or the vendored
``aslscan-Cargo.lock`` in ``aslprep/tests/aslscan_fixtures.py``, regenerate the spec file, and
rerun the integration tier; record any moved bound in the same pull request. Phantoms are
pluggable (``PHANTOMS`` in the spec): the TemplateFlow-based ``tfmni`` phantoms are built
in place, and a future ``subject:<id>`` phantom would download a real subject's CBF, ATT,
aBAT, aBV, T1, T2 and T2* maps; ``sanitize_maps`` already makes such maps satisfy aslscan's
phantom contract.

If a cached CircleCI fixture entry is bad, bump ``CACHE_EPOCH`` and regenerate the spec file;
CircleCI caches cannot be overwritten under the same key.
