"""Interfaces for partial volume correction."""

import nibabel as nb
import nitransforms as nt
import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep import config
from aslprep.utils.pvc import (
    PVC_STATUS,
    build_pv_maps,
    calibrate_csf,
    calibrate_pvc,
    get_slice_axis,
    kernel_size_to_axes,
    local_regression_csf,
    local_regression_pvc,
    normalized_resample_ratio,
    resolve_effective_resolution,
    summarize_csf_qc,
)


class _PartialVolumeMapsInputSpec(BaseInterfaceInputSpec):
    gm_tpm = File(exists=True, mandatory=True, desc='Gray matter probability map (anatomical).')
    wm_tpm = File(exists=True, mandatory=True, desc='White matter probability map (anatomical).')
    csf_tpm = File(exists=True, mandatory=True, desc='CSF probability map (anatomical).')
    reference = File(
        exists=True,
        mandatory=True,
        desc='An image on the ASL reference grid, such as the ASL brain mask.',
    )
    aslref2anat_xfm = File(
        exists=True,
        mandatory=True,
        desc='Affine transform from the ASL reference to the anatomical image (ITK format).',
    )
    metadata = traits.Dict(mandatory=True, desc='ASL metadata.')
    effective_resolution = traits.Either(
        traits.Enum('nominal', 'auto'),
        traits.List(traits.Float, minlen=3, maxlen=3),
        default='nominal',
        usedefault=True,
        desc=(
            "'nominal', 'auto', or the total response FWHM in mm, "
            'ordered (in-plane, in-plane, slice).'
        ),
    )
    max_points = traits.Int(
        8_000_000,
        usedefault=True,
        nohash=True,
        desc='Supersampled points processed at once. Bounds the memory used.',
    )


class _PartialVolumeMapsOutputSpec(TraitedSpec):
    gm_pv = File(exists=True, desc='Gray matter PV map on the ASL reference grid.')
    wm_pv = File(exists=True, desc='White matter PV map on the ASL reference grid.')
    csf_pv = File(exists=True, desc='CSF PV map on the ASL reference grid.')
    pv_metadata = traits.Dict(desc='How the PV maps were built, for their sidecars.')


class PartialVolumeMaps(SimpleInterface):
    """Build tissue PV maps on the ASL grid that model its effective resolution.

    See :func:`~aslprep.utils.pvc.build_pv_maps` and
    :func:`~aslprep.utils.pvc.resolve_effective_resolution`.
    """

    input_spec = _PartialVolumeMapsInputSpec
    output_spec = _PartialVolumeMapsOutputSpec

    def _run_interface(self, runtime):
        ref_img = nb.load(self.inputs.reference)
        ref_affine = ref_img.affine
        voxel_size = tuple(float(v) for v in ref_img.header.get_zooms()[:3])
        slice_axis = get_slice_axis(self.inputs.metadata)

        resolution = resolve_effective_resolution(
            self.inputs.effective_resolution,
            self.inputs.metadata,
            voxel_size,
            slice_axis,
        )
        for warning in resolution['warnings']:
            config.loggers.interface.warning(warning)

        # The aslref-to-anat transform maps anatomical points to aslref points (it pulls the
        # anatomical image onto the ASL grid), so its inverse takes aslref points to anat.
        xfm = nt.linear.load(self.inputs.aslref2anat_xfm)
        matrix = np.asarray(xfm.matrix)
        if matrix.ndim == 3:
            if matrix.shape[0] != 1:
                raise ValueError('Expected a single affine transform from aslref to anat.')
            matrix = matrix[0]
        ref2anat = np.linalg.inv(matrix)

        tpm_imgs = [
            nb.load(f) for f in (self.inputs.gm_tpm, self.inputs.wm_tpm, self.inputs.csf_tpm)
        ]
        if not all(np.allclose(img.affine, tpm_imgs[0].affine) for img in tpm_imgs[1:]):
            raise ValueError('The tissue probability maps must share a grid.')

        pv_maps, info = build_pv_maps(
            [img.get_fdata(dtype=np.float32) for img in tpm_imgs],
            tpm_imgs[0].affine,
            ref_img.shape[:3],
            ref_affine,
            ref2anat,
            fwhm=resolution['fwhm'],
            max_points=self.inputs.max_points,
            slice_axis=slice_axis,
        )

        for name, data in zip(('gm_pv', 'wm_pv', 'csf_pv'), pv_maps, strict=True):
            out_img = nb.Nifti1Image(data.astype(np.float32), ref_affine, ref_img.header)
            out_img.set_data_dtype(np.float32)
            self._results[name] = fname_presuffix(
                self.inputs.reference,
                suffix=f'_{name}',
                newpath=runtime.cwd,
                use_ext=True,
            )
            out_img.to_filename(self._results[name])

        self._results['pv_metadata'] = {
            'EffectiveResolution': [float(v) for v in resolution['fwhm']],
            'EffectiveResolutionSource': resolution['source'],
            'EffectiveResolutionFactors': resolution['factors'],
            'EffectiveResolutionMetadata': resolution['metadata_used'],
            'SliceEncodingAxis': 'ijk'[slice_axis],
            'Supersampling': info['subdivision'],
            'GaussianSigma': info['sigma_mm'],
            'Description': (
                'Tissue fraction on the ASL reference grid: the tissue probability map, '
                'linearly interpolated, averaged over each voxel and convolved with a Gaussian '
                'so that the total response has the effective resolution (FWHM, mm, per voxel '
                'axis). Blurring from motion correction and distortion correction resampling '
                'is not modelled.'
            ),
        }
        return runtime


class _LocalRegressionPVCInputSpec(BaseInterfaceInputSpec):
    cbf = File(exists=True, mandatory=True, desc='Mean CBF map, quantified by ComputeCBF.')
    m0_calib = File(exists=True, mandatory=True, desc='The M0 image CBF was divided by.')
    quant_valid = File(exists=True, mandatory=True, desc='Where CBF could be quantified.')
    mask = File(exists=True, mandatory=True, desc='Brain mask.')
    gm_pv = File(exists=True, mandatory=True, desc='Gray matter PV map.')
    wm_pv = File(exists=True, mandatory=True, desc='White matter PV map.')
    csf_pv = File(
        exists=True,
        mandatory=True,
        desc='CSF PV map, for the CSF quality-control fit (the GM/WM fit assumes CSF is 0).',
    )
    metadata = traits.Dict(mandatory=True, desc='ASL metadata, for the slice-encoding axis.')
    kernel_size = traits.List(
        traits.Int,
        [5, 5, 1],
        minlen=3,
        maxlen=3,
        usedefault=True,
        desc='Kernel size in voxels, ordered (in-plane, in-plane, slice).',
    )
    tau = traits.Float(0.5, usedefault=True, desc='Minimum tissue support in the kernel.')
    kappa_max = traits.Float(
        100.0, usedefault=True, desc='Maximum condition number of the normalized Gram matrix.'
    )
    c_max = traits.Float(
        0.05, usedefault=True, desc='Maximum contamination ratio for single-tissue fits.'
    )
    min_pv = traits.Float(
        0.1, usedefault=True, desc="A tissue's estimate is written only where its PV is this high."
    )


class _LocalRegressionPVCOutputSpec(TraitedSpec):
    gm_cbf = File(exists=True, desc='Gray matter CBF.')
    wm_cbf = File(exists=True, desc='White matter CBF.')
    gm_se = File(exists=True, desc='Standard error of gray matter CBF.')
    wm_se = File(exists=True, desc='Standard error of white matter CBF.')
    gm_valid = File(exists=True, desc='Where gray matter CBF was estimated.')
    wm_valid = File(exists=True, desc='Where white matter CBF was estimated.')
    status = File(exists=True, desc='Fit status of each voxel.')
    condition = File(exists=True, desc='Condition number of the normalized Gram matrix.')
    gm_noise = File(exists=True, desc='Noise amplification of the gray matter estimate.')
    wm_noise = File(exists=True, desc='Noise amplification of the white matter estimate.')
    rmse = File(exists=True, desc='Residual root mean square error of each fit.')
    csf_cbf = File(exists=True, desc='CSF "CBF" from the three-tissue QC fit (expected near 0).')
    csf_se = File(exists=True, desc='Standard error of the CSF estimate.')
    csf_status = File(exists=True, desc='Status of the three-tissue CSF QC fit.')
    csf_noise = File(exists=True, desc='Noise amplification of the CSF estimate.')
    summary = traits.Dict(desc='Fraction of brain voxels with each status.')
    qc = traits.Dict(desc='Scalar quality-control measures, for the CBF QC file.')
    pvc_metadata = traits.Dict(desc='Settings, for the sidecars.')


class LocalRegressionPVC(SimpleInterface):
    """Partial volume correction by local linear regression.

    The M0-free signal, CBF times the M0 it was quantified with, is regressed on gray and
    white matter PV maps within a kernel around each voxel
    :footcite:p:`asllani2008regression`, and the tissue amplitudes are calibrated with the
    centre voxel's M0. See :func:`~aslprep.utils.pvc.local_regression_pvc`.

    References
    ----------
    .. footbibliography::
    """

    input_spec = _LocalRegressionPVCInputSpec
    output_spec = _LocalRegressionPVCOutputSpec

    def _run_interface(self, runtime):
        cbf_img = nb.load(self.inputs.cbf)
        cbf = cbf_img.get_fdata()
        m0 = nb.load(self.inputs.m0_calib).get_fdata()
        quant_valid = nb.load(self.inputs.quant_valid).get_fdata() > 0
        mask = nb.load(self.inputs.mask).get_fdata() > 0
        pv_gm = nb.load(self.inputs.gm_pv).get_fdata()
        pv_wm = nb.load(self.inputs.wm_pv).get_fdata()

        slice_axis = get_slice_axis(self.inputs.metadata)
        kernel = kernel_size_to_axes(self.inputs.kernel_size, slice_axis)

        # The M0-free signal: q * delta-M
        signal = np.where(quant_valid, cbf * m0, np.nan)
        fit = local_regression_pvc(
            signal,
            pv_gm,
            pv_wm,
            mask=mask,
            valid=quant_valid,
            kernel_size=kernel,
            tau=self.inputs.tau,
            kappa_max=self.inputs.kappa_max,
            c_max=self.inputs.c_max,
        )
        calibrated = calibrate_pvc(fit, m0, pv_gm, pv_wm, min_pv=self.inputs.min_pv)

        # Quality control: a separate fit with CSF as a third tissue. CSF has no perfusion,
        # so its estimate should be near zero. The GM and WM maps above are unaffected.
        pv_csf = nb.load(self.inputs.csf_pv).get_fdata()
        csf_fit = local_regression_csf(
            signal,
            pv_gm,
            pv_wm,
            pv_csf,
            mask=mask,
            valid=quant_valid,
            kernel_size=kernel,
            tau=self.inputs.tau,
            kappa_max=self.inputs.kappa_max,
            c_max=self.inputs.c_max,
        )
        csf = calibrate_csf(csf_fit, m0, pv_csf, min_pv=self.inputs.min_pv)
        csf_summary = summarize_csf_qc(csf['cbf_csf'], csf['valid_csf'], pv_csf)

        outputs = {
            'gm_cbf': (calibrated['cbf_gm'], np.float32),
            'wm_cbf': (calibrated['cbf_wm'], np.float32),
            'gm_se': (calibrated['se_gm'], np.float32),
            'wm_se': (calibrated['se_wm'], np.float32),
            'gm_valid': (calibrated['valid_gm'], np.uint8),
            'wm_valid': (calibrated['valid_wm'], np.uint8),
            'status': (fit['status'], np.uint8),
            'condition': (fit['condition'], np.float32),
            'gm_noise': (fit['noise_gm'], np.float32),
            'wm_noise': (fit['noise_wm'], np.float32),
            'rmse': (fit['rmse'], np.float32),
            'csf_cbf': (csf['cbf_csf'], np.float32),
            'csf_se': (csf['se_csf'], np.float32),
            'csf_status': (csf_fit['status'], np.uint8),
            'csf_noise': (np.where(csf['valid_csf'], csf_fit['noise_csf'], 0.0), np.float32),
        }
        for name, (data, dtype) in outputs.items():
            out_img = nb.Nifti1Image(data.astype(dtype), cbf_img.affine, cbf_img.header)
            out_img.set_data_dtype(dtype)
            self._results[name] = fname_presuffix(
                self.inputs.cbf,
                suffix=f'_lr{name}',
                newpath=runtime.cwd,
            )
            out_img.to_filename(self._results[name])

        n_mask = max(int(mask.sum()), 1)
        self._results['summary'] = {
            PVC_STATUS[code]: float(np.sum(fit['status'][mask] == code) / n_mask)
            for code in (1, 2, 3, 4)
        }
        self._results['qc'] = {
            'lr_csf_cbf_median': csf_summary['median'],
            'lr_csf_cbf_iqr': csf_summary['iqr'],
            'lr_csf_n_voxels': csf_summary['n'],
            'lr_fraction_two_tissue': self._results['summary'][PVC_STATUS[1]],
        }
        self._results['pvc_metadata'] = {
            'KernelSize': list(self.inputs.kernel_size),
            'KernelAxes': ['in-plane', 'in-plane', 'slice'],
            'MinimumTissueSupport': self.inputs.tau,
            'MaximumConditionNumber': self.inputs.kappa_max,
            'MaximumContamination': self.inputs.c_max,
            'MinimumPartialVolume': self.inputs.min_pv,
            'Calibration': (
                "Each tissue's signal amplitude is divided by the centre voxel's M0, the same "
                'M0 the uncorrected CBF was quantified with.'
            ),
        }
        return runtime


class _DivideByWeightInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='The warped (value * validity) image.')
    weight_file = File(exists=True, mandatory=True, desc='The warped validity image.')
    threshold = traits.Float(
        0.5, usedefault=True, desc='Voxels with less warped validity than this are set to zero.'
    )


class _DivideByWeightOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='The validity-normalized image.')


class DivideByWeight(SimpleInterface):
    """Finish a validity-normalized resampling by dividing by the resampled validity.

    Zeros that stand for "not estimated" would otherwise be blended into the resampled values.
    See :func:`~aslprep.utils.pvc.normalized_resample_ratio`.
    """

    input_spec = _DivideByWeightInputSpec
    output_spec = _DivideByWeightOutputSpec

    def _run_interface(self, runtime):
        img = nb.load(self.inputs.in_file)
        weight = nb.load(self.inputs.weight_file).get_fdata()
        data = normalized_resample_ratio(img.get_fdata(), weight, self.inputs.threshold)
        out_img = nb.Nifti1Image(data.astype(np.float32), img.affine, img.header)
        out_img.set_data_dtype(np.float32)
        self._results['out_file'] = fname_presuffix(
            self.inputs.in_file, suffix='_normalized', newpath=runtime.cwd
        )
        out_img.to_filename(self._results['out_file'])
        return runtime
