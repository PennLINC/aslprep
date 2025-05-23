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

   - This value (especially for the M0 scan) may be useful for calculating CBF.
     It's not used in ASLPrep, but other tools may use it.

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

The M0 scan is not necessarily used for distortion correction.
It will be used for CBF calculations, and may be used for other purposes within ASLPrep
(e.g., as a high-contrast image for coregistration).


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
- Potentially others.
