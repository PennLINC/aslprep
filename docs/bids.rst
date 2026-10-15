.. include:: links.rst

###################
Preparing BIDS Data
###################

Here we discuss a few different quirks in BIDS that might affect how you prepare your ASL data.


***********************
Mixed-type 4D ASL files
***********************

Unlike fMRI data, the 4D ASL file can contain multiple types of volumes, including M0, control,
label, noise, and derived volumes.
As such, the time axis is not as meaningful as it is for fMRI data.

The aslcontext.tsv file is used to describe the volume types in the ASL file.

The NIfTI file format only allows for a single time-step unit,
even though ASL files may have different time-steps for different volumes.
This has necessitated metadata in the ASL NIfTI files' sidecar JSON files to describe the time-steps
for each volume; namely RepetitionTimePreparation and RepetitionTimeExcitation.


***********************
Repetition Time Nuances
***********************

.. warning::

   It's important to talk to your physicist to determine the correct values for these fields.
   The values produced by the DICOM converter are not always correct- probably because the DICOM
   does not have the right information.

There are three relevant repetition time fields in ASL data:

1. RepetitionTime: This can only be a single value. When the ASL BEP was first developed,
   it was decided that the RepetitionTime field would not be changed to allow arrays of values,
   as is needed for ASL data. Instead, they created the RepetitionTimePreparation field, which
   can be either a single value or an array of values (one value per volume in the ASL file).
2. RepetitionTimePreparation: This is effectively the same as RepetitionTime, but can be an array.

   - ASLPrep uses this value for the M0 volumes (or the control volumes, when there is no M0)
     to scale the M0 data when the repetition time is shorter than 5 seconds.
     It must be greater than zero.
     ASLPrep raises an error otherwise, because CBF cannot be calculated.

3. RepetitionTimeExcitation: The "true" time between the beginning of the first excitation pulse
   and the last one.


**********
M0 Scaling
**********

In many ASL sequences, the M0 volumes need to be scaled before they can be used to calculate CBF.
This is one value that many dataset curators will not know offhand.
Unfortunately, there is no standard way to find this value.
Some sequences will include the scaling factor in a private DICOM field,
or it may be available in the protocol's summary file (see below for an example).

BIDS expects the M0 volumes to be scaled as part of the BIDS curation process,
but ASLPrep has a ``--m0-scale`` parameter that users can use if they haven't already scaled
the M0 volumes.

.. tip::

   If you are working with a BIDS dataset where it is not clear if the M0 volumes have been scaled,
   and, if they have not, it's not clear how much they need to be scaled, you can compute CBF based
   on M0 scale = 1 and then try to figure out the actual value if the values look out of range.


*********************************************
Separate M0 scans and the "IntendedFor" field
*********************************************

When you have a separate M0 scan, you need to explicitly link it to the ASL series with which it will be used.
The M0 scan can be linked to multiple ASL series, or each ASL series can have its own M0 scan.

The ``IntendedFor`` field is used to link the M0 scan to the ASL series,
but it is not used in the same way as other ``IntendedFor`` cases (e.g., to link field maps to images).

Each entry in ``IntendedFor`` may be a BIDS URI, relative to the dataset root,
as recommended since BIDS 1.7,
or a path relative to the subject directory, which BIDS has deprecated but still allows.
For example, in ``sub-01/ses-1/perf/sub-01_ses-1_m0scan.json``:

.. code-block:: json

   {
       "IntendedFor": ["bids::sub-01/ses-1/perf/sub-01_ses-1_asl.nii.gz"]
   }

is equivalent to ``"IntendedFor": ["ses-1/perf/sub-01_ses-1_asl.nii.gz"]``.

The M0 scan is not necessarily used for distortion correction.
It will be used for CBF calculations, and may be used for other purposes within ASLPrep
(e.g., as a high-contrast image for coregistration).


.. _haematocrit:

**********************************
Haematocrit and the T1 of blood
**********************************

CBF depends directly on the T1 of arterial blood (:math:`T_{1,blood}`),
which in turn depends on the haematocrit: the fraction of the blood's volume taken up by
red blood cells.
By default, *ASLPrep* uses a fixed :math:`T_{1,blood}` for the scanner's field strength
(1.65 s at 3 T :footcite:p:`alsop_recommended_2015`), which assumes a typical haematocrit.
In people whose haematocrit is far from typical, such as people with anaemia or sickle cell
disease, that assumption can bias CBF by tens of percent :footcite:p:`hales2016general`.

*ASLPrep* will automatically use a measured haematocrit when the dataset includes one.
For each ASL file, these sources are checked in order, and the first usable value is used:

1. **ASL sidecar:** a ``Haematocrit`` field in the ASL file's JSON sidecar
   (including sidecars inherited from higher levels of the dataset).
2. **Scans file:** ``<bids-root>/sub-<label>/[ses-<label>/]sub-<label>[_ses-<label>]_scans.tsv``,
   in the row for the ASL file.
3. **Sessions file:** ``<bids-root>/sub-<label>/sub-<label>_sessions.tsv``,
   in the row for the ASL file's session.
4. **Participants file:** ``<bids-root>/participants.tsv``, in the row for the participant.

In the TSV files, the value goes in a ``haematocrit`` column.
The American spelling (``Hematocrit`` or ``hematocrit``) is also accepted, in any of the sources.
A missing value (``n/a``) moves on to the next source.

The value should be a fraction between 0 and 1, as BIDS defines it.
Values between 1 and 100 are read as percentages, and other values are ignored with a warning.
We recommend adding a JSON sidecar for the TSV file that describes the column, for example::

   {
       "haematocrit": {
           "Description": "Haematocrit measured from a blood sample on the day of the scan",
           "Units": "fraction"
       }
   }

.. note::

   BIDS does not yet define a haematocrit field for ASL data.
   ``Haematocrit`` is the field BIDS defines for PET blood data (``*_blood.json``),
   so *ASLPrep* uses the same name, spelling and definition.
   Extra columns in the participants, sessions and scans files are allowed in BIDS.

When a haematocrit is found, *ASLPrep* scales its usual :math:`T_{1,blood}` for that
haematocrit, using the model of :footcite:t:`hales2016general`.
The model accounts for the haematocrit, the field strength and the oxygen saturation of the
blood (assumed to be 0.97, for arterial blood):

.. math::

   T_{1,blood} = T_{1,blood}^{default} \cdot
   \frac{T_{1,blood}^{Hales}(\mathrm{Hct})}{T_{1,blood}^{Hales}(0.42)}

The default values are taken to correspond to a typical haematocrit of 0.42
:footcite:p:`lu2004determining`, so a participant with that haematocrit gets exactly the
default, and CBF stays comparable between participants with and without a measured haematocrit.
The published model is not used on its own, because at a typical haematocrit it gives values
7-21% longer than the defaults, which would shift CBF for every participant with a measurement.

.. list-table:: :math:`T_{1,blood}` (s) for a few haematocrit values
   :header-rows: 1

   * - Haematocrit
     - 1.5 T
     - 3 T
     - 7 T
   * - 0.20
     - 1.596
     - 1.948
     - 2.447
   * - 0.30
     - 1.477
     - 1.804
     - 2.274
   * - 0.42 (default)
     - 1.350
     - 1.650
     - 2.087
   * - 0.47
     - 1.301
     - 1.591
     - 2.015

For single-delay PCASL at 3 T (1.8 s labeling duration and post-labeling delay),
a haematocrit of 0.20 instead of a typical one lowers CBF by about 21%.
The same :math:`T_{1,blood}` is used for BASIL, and the methods boilerplate reports it.

The model is valid between 1.5 and 7 T.
At other field strengths, the haematocrit is ignored with a warning,
and :math:`T_{1,blood}` is set from the field strength as usual.
Without a haematocrit, *ASLPrep* behaves as before.

.. important::

   *ASLPrep* does not estimate the haematocrit from the participant's age or sex.
   Population-average haematocrit differs by age and sex,
   but using an average in place of a measurement would change CBF for every dataset with
   ``age`` and ``sex`` columns, without improving accuracy for the participants whose
   haematocrit is atypical.
   If you want to use population-average values, add them to the ``haematocrit`` column yourself.

.. footbibliography::


*********************
Other metadata fields
*********************

One of the most common issues with working with ASL data is that many of the required metadata fields
are not encoded in standard DICOM fields, so they cannot be inferred from the DICOMs alone
(at least not by a layperson).
We recommend checking with your scanner physicist, when possible, to determine the correct values
for these fields.

Ask your physicist for the correct values for the following fields:

- If the sequence is PCASL or PASL. This is typically in the series description or protocol name, but it can be wrong.
- ``BackgroundSuppression``
- ``LabelingEfficiency`` (or rely on ASLPrep's heuristic)
- ``PostLabelingDelay`` (should be in protocol summary file and private DICOM fields, but may be difficult to discover)
- ``LabelingDuration`` (should be in protocol summary file and private DICOM fields, but may be difficult to discover)
- ``BackgroundSuppressionNumberPulses``
- ``VascularCrushing``
- ``LookLocker``
- The order of volumes in the ASL file.
  The aslcontext.tsv is not automatically generated and is not something that can be inferred
  from the DICOMs (at least not easily).

The list above is not comprehensive.
There may be other fields or elements that need to be manually added to the dataset before it is usable.
