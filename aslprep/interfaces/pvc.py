"""Interfaces for partial volume correction."""

import nibabel as nb
import nitransforms as nt
import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep import config
from aslprep.utils.pvc import (
    CSF_QC_DEFAULT_ATT,
    LR_ABV_MAX,
    LR_CBF_MAX,
    PVC_STATUS,
    build_pv_maps,
    calibrate_csf,
    calibrate_pvc,
    csf_kinetic_projection,
    fit_tissue_curves,
    get_slice_axis,
    gkm_predict,
    kernel_size_to_axes,
    local_regression_csf,
    local_regression_pvc,
    normalized_resample_ratio,
    resolve_effective_resolution,
    summarize_csf_qc,
)

#: PV below which a tissue's T1 on the ASL grid is not defined (and the scalar T1 is used).
MIN_T1_WEIGHT = 1e-6


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
    gm_r1 = File(
        exists=True,
        desc=(
            'Gray matter R1 (1/s) in each anatomical voxel. If given with wm_r1, the tissue T1s '
            'are averaged onto the ASL reference grid with the PV maps.'
        ),
    )
    wm_r1 = File(exists=True, desc='White matter R1 (1/s) in each anatomical voxel.')


class _PartialVolumeMapsOutputSpec(TraitedSpec):
    gm_pv = File(exists=True, desc='Gray matter PV map on the ASL reference grid.')
    wm_pv = File(exists=True, desc='White matter PV map on the ASL reference grid.')
    csf_pv = File(exists=True, desc='CSF PV map on the ASL reference grid.')
    pv_metadata = traits.Dict(desc='How the PV maps were built, for their sidecars.')
    gm_t1 = File(
        exists=True,
        desc=(
            'Gray matter T1 (s) on the ASL reference grid: 1 over the mean gray matter R1, '
            'weighted as the PV map weights the tissue. Zero where there is no gray matter.'
        ),
    )
    wm_t1 = File(exists=True, desc='White matter T1 (s) on the ASL reference grid, likewise.')


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
        # The same voxel size build_pv_maps integrates over (the affine's column norms), so
        # 'nominal' always means no blur beyond the voxel, even if the header zooms disagree.
        voxel_size = tuple(float(v) for v in np.sqrt((ref_affine[:3, :3] ** 2).sum(axis=0)))
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

        maps = [img.get_fdata(dtype=np.float32) for img in tpm_imgs]
        with_t1 = isdefined(self.inputs.gm_r1) and isdefined(self.inputs.wm_r1)
        if with_t1:
            # Each tissue's R1, weighted by its probability, is averaged with the same response
            # as the PV map, so that the mean R1 weights the anatomical voxels as the PV map does.
            for tpm, r1_file in zip(maps[:2], (self.inputs.gm_r1, self.inputs.wm_r1), strict=True):
                r1_img = nb.load(r1_file)
                if r1_img.shape[:3] != tpm.shape or not np.allclose(
                    r1_img.affine,
                    tpm_imgs[0].affine,
                    atol=1e-3,
                ):
                    raise ValueError(f'{r1_file} is not on the grid of the tissue maps.')

                maps.append(np.clip(tpm, 0, 1) * r1_img.get_fdata(dtype=np.float32))

        pv_maps, info = build_pv_maps(
            maps,
            tpm_imgs[0].affine,
            ref_img.shape[:3],
            ref_affine,
            ref2anat,
            fwhm=resolution['fwhm'],
            max_points=self.inputs.max_points,
            slice_axis=slice_axis,
            clip=[True, True, True] + [False] * (len(maps) - 3),
        )

        if with_t1:
            for name, pv, weighted in (
                ('gm_t1', pv_maps[0], pv_maps[3]),
                ('wm_t1', pv_maps[1], pv_maps[4]),
            ):
                with np.errstate(divide='ignore', invalid='ignore'):
                    t1 = np.where((pv > MIN_T1_WEIGHT) & (weighted > 0), pv / weighted, 0.0)
                out_img = nb.Nifti1Image(t1.astype(np.float32), ref_affine, ref_img.header)
                out_img.set_data_dtype(np.float32)
                self._results[name] = fname_presuffix(
                    self.inputs.reference,
                    suffix=f'_{name}',
                    newpath=runtime.cwd,
                    use_ext=True,
                )
                out_img.to_filename(self._results[name])

        for name, data in zip(('gm_pv', 'wm_pv', 'csf_pv'), pv_maps[:3], strict=True):
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


class _MultiDelayLocalRegressionPVCInputSpec(BaseInterfaceInputSpec):
    deltam = File(exists=True, mandatory=True, desc='Delta-M series, as ExtractCBF writes it.')
    metadata = traits.Dict(
        mandatory=True,
        desc='ASL metadata of the delta-M series (after volume selection), with its delays.',
    )
    calibration = traits.Dict(
        desc=(
            'Calibration specification, from aslprep.utils.calibration.resolve_calibration. '
            'If not defined, it is resolved from metadata.'
        ),
    )
    m0_calib = File(exists=True, mandatory=True, desc='The M0 image the fit is calibrated with.')
    quant_valid = File(exists=True, mandatory=True, desc='Where the data can be fit.')
    mask = File(exists=True, mandatory=True, desc='Brain mask.')
    gm_pv = File(exists=True, mandatory=True, desc='Gray matter PV map.')
    wm_pv = File(exists=True, mandatory=True, desc='White matter PV map.')
    csf_pv = File(exists=True, mandatory=True, desc='CSF PV map, for the CSF QC fit.')
    mean_cbf = File(
        exists=True,
        mandatory=True,
        desc='Uncorrected CBF: where it is positive, the uncorrected ATT was estimated.',
    )
    att = File(
        exists=True,
        mandatory=True,
        desc='Uncorrected ATT, which sets the shape of the CSF quality-control basis.',
    )
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
    num_threads = traits.Int(1, usedefault=True, nohash=True, desc='Threads for the kinetic fits.')
    use_gpu = traits.Bool(False, usedefault=True, desc='Run the kinetic fits on a GPU, with CuPy.')
    gm_t1 = File(
        exists=True,
        desc=(
            'Gray matter T1 (s) in each voxel, for the tissue-T1 kinetic model '
            "(--t1map-method local-regression). Where it is zero, and without it, the calibration's "
            'GrayMatterT1 is used.'
        ),
    )
    wm_t1 = File(exists=True, desc='White matter T1 (s) in each voxel, likewise.')


class _MultiDelayLocalRegressionPVCOutputSpec(TraitedSpec):
    gm_cbf = File(exists=True, desc='Gray matter CBF.')
    wm_cbf = File(exists=True, desc='White matter CBF.')
    gm_valid = File(exists=True, desc='Where gray matter CBF was estimated.')
    wm_valid = File(exists=True, desc='Where white matter CBF was estimated.')
    gm_att = File(exists=True, desc='Gray matter arterial transit time.')
    wm_att = File(exists=True, desc='White matter arterial transit time.')
    gm_att_valid = File(exists=True, desc='Where gray matter ATT was estimated.')
    wm_att_valid = File(exists=True, desc='Where white matter ATT was estimated.')
    gm_att_lower = File(exists=True, desc='Shortest gray matter ATT that fits as well.')
    gm_att_upper = File(exists=True, desc='Longest gray matter ATT that fits as well.')
    wm_att_lower = File(exists=True, desc='Shortest white matter ATT that fits as well.')
    wm_att_upper = File(exists=True, desc='Longest white matter ATT that fits as well.')
    gm_kinstatus = File(exists=True, desc='Kinetic fit status of the gray matter curve.')
    wm_kinstatus = File(exists=True, desc='Kinetic fit status of the white matter curve.')
    gm_kinres = File(exists=True, desc='Relative kinetic residual of the gray matter curve.')
    wm_kinres = File(exists=True, desc='Relative kinetic residual of the white matter curve.')
    status = File(exists=True, desc='Spatial regression status of each voxel.')
    condition = File(exists=True, desc='Condition number of the normalized Gram matrix.')
    gm_noise = File(exists=True, desc='Noise amplification of the gray matter amplitudes.')
    wm_noise = File(exists=True, desc='Noise amplification of the white matter amplitudes.')
    rmse = File(exists=True, desc='Pooled spatial residual of each regression, in delta-M units.')
    csf_cbf = File(exists=True, desc='Signed apparent CSF CBF (QC, expected near 0).')
    csf_status = File(exists=True, desc='Status of the three-tissue CSF QC fit.')
    csf_noise = File(exists=True, desc='Noise amplification of the CSF amplitudes.')
    csf_kinres = File(exists=True, desc='Relative residual of the CSF curve against its template.')
    summary = traits.Dict(desc='Fraction of brain voxels with each status.')
    qc = traits.Dict(desc='Scalar quality-control measures, for the CBF QC file.')
    pvc_metadata = traits.Dict(desc='Settings, for the sidecars.')


class MultiDelayLocalRegressionPVC(SimpleInterface):
    """Partial volume correction of multi-delay data by local regression (spec Section 5.7).

    Each delay volume's delta-M is unmixed into gray and white matter amplitudes with one
    spatial design (:func:`~aslprep.utils.pvc.local_regression_pvc`), and each tissue's curve is
    then fit with the general kinetic model, calibrated by the centre voxel's M0
    (:func:`~aslprep.utils.pvc.fit_tissue_curves`) :footcite:p:`asllani2008regression`.

    References
    ----------
    .. footbibliography::
    """

    input_spec = _MultiDelayLocalRegressionPVCInputSpec
    output_spec = _MultiDelayLocalRegressionPVCOutputSpec

    def _run_interface(self, runtime):
        from nilearn import maskers

        from aslprep.utils.asl import pcasl_or_pasl
        from aslprep.utils.calibration import (
            calibration_value,
            m0_partition_coefficient,
            resolve_calibration,
            tissue_t1_kinetics,
            voxelwise_tissue_t1prime,
        )
        from aslprep.utils.cbf import (
            get_array_module,
            multipld_bolus,
            voxel_delays,
        )

        metadata = self.inputs.metadata
        deltam_img = nb.load(self.inputs.deltam)
        deltam = deltam_img.get_fdata()
        if deltam.ndim != 4:
            raise ValueError(f'Expected a 4D delta-M series, got shape {deltam.shape}.')
        shape = deltam.shape[:3]
        n_vol = deltam.shape[3]
        affine = deltam_img.affine

        m0 = nb.load(self.inputs.m0_calib).get_fdata()
        quant_valid = nb.load(self.inputs.quant_valid).get_fdata() > 0
        mask = nb.load(self.inputs.mask).get_fdata() > 0
        pv_gm = nb.load(self.inputs.gm_pv).get_fdata()
        pv_wm = nb.load(self.inputs.wm_pv).get_fdata()
        pv_csf = nb.load(self.inputs.csf_pv).get_fdata()

        is_casl = pcasl_or_pasl(metadata=metadata)
        if isdefined(self.inputs.calibration):
            calibration = self.inputs.calibration
        else:
            calibration = resolve_calibration(metadata)
        tau, ti1 = multipld_bolus(metadata, is_casl)
        # The same parameters as the uncorrected fit (ComputeCBF), which made m0_calib
        kinetics = {
            'labeleff': calibration_value(calibration, 'EffectiveLabelingEfficiency'),
            't1blood': calibration_value(calibration, 'BloodT1'),
            'partition_coefficient': m0_partition_coefficient(calibration),
            'is_casl': is_casl,
            'tau': tau,
            'ti1': ti1,
        }

        # Delays of every voxel and volume, in the same way ComputeCBF builds them
        all_voxels = nb.Nifti1Image(np.ones(shape, np.uint8), affine)
        masker = maskers.NiftiMasker(mask_img=all_voxels).fit()
        plds, _ = voxel_delays(metadata, deltam_img, masker, int(np.prod(shape)), True)
        plds = masker.inverse_transform(plds.T).get_fdata()
        if plds.ndim == 3:
            plds = plds[..., None]
        if plds.shape[3] != n_vol:
            raise ValueError(
                f'Number of delays ({plds.shape[3]}) does not match delta-M volumes ({n_vol}).'
            )

        slice_axis = get_slice_axis(metadata)
        kernel = kernel_size_to_axes(self.inputs.kernel_size, slice_axis)
        lr_kwargs = {
            'kernel_size': kernel,
            'tau': self.inputs.tau,
            'kappa_max': self.inputs.kappa_max,
            'c_max': self.inputs.c_max,
        }

        # Step 1: one spatial design for every volume
        signal = np.where(quant_valid[..., None], deltam, np.nan)
        fit = local_regression_pvc(signal, pv_gm, pv_wm, mask=mask, valid=quant_valid, **lr_kwargs)
        calibratable = np.isfinite(m0) & (m0 > 0)
        spatial = {
            'gm': np.isin(fit['status'], (1, 2)) & (pv_gm >= self.inputs.min_pv) & calibratable,
            'wm': np.isin(fit['status'], (1, 3)) & (pv_wm >= self.inputs.min_pv) & calibratable,
        }

        # Step 2: a kinetic fit of each tissue's curve, calibrated by the center voxel's M0.
        # With the tissue-T1 model, each curve decays with its own tissue's apparent T1.
        tissue_t1 = tissue_t1_kinetics(calibration)
        t1prime = None if tissue_t1 is None else tissue_t1['t1prime']
        if tissue_t1 is not None and isdefined(self.inputs.gm_t1) and isdefined(self.inputs.wm_t1):
            tissue_t1s = []
            for t1_file in (self.inputs.gm_t1, self.inputs.wm_t1):
                t1_img = nb.load(t1_file)
                if t1_img.shape[:3] != shape or not np.allclose(t1_img.affine, affine, atol=1e-3):
                    raise ValueError(f'{t1_file} does not match the delta-M grid ({shape}).')

                tissue_t1s.append(t1_img.get_fdata().reshape(-1))

            # One apparent T1 per voxel and tissue, on the image grid
            t1prime = voxelwise_tissue_t1prime(calibration, np.column_stack(tissue_t1s))
            t1prime = t1prime.reshape(*shape, 2)

        xp = get_array_module(self.inputs.use_gpu)
        maps = {}
        for i_tissue, tissue in enumerate(('gm', 'wm')):
            sel = spatial[tissue]
            tissue_t1prime = None
            if t1prime is not None:
                # Each curve is of one tissue, so it decays with that tissue's T1' alone
                tissue_t1prime = (
                    t1prime[i_tissue] if t1prime.ndim == 1 else t1prime[sel][:, i_tissue, None]
                )

            result = fit_tissue_curves(
                fit[f'beta_{tissue}'][sel],
                m0[sel],
                plds[sel],
                xp=xp,
                n_jobs=self.inputs.num_threads,
                t1prime=tissue_t1prime,
                **kinetics,
            )
            status = np.zeros(shape, dtype=np.uint8)
            status[sel] = result['status']
            cbf_valid = np.isin(status, (1, 2))
            att_valid = status == 1
            for key in ('cbf', 'att', 'att_lower', 'att_upper', 'kinres'):
                arr = np.zeros(shape)
                arr[sel] = result[key]
                maps[f'{tissue}_{key}'] = arr
            maps[f'{tissue}_cbf'] = np.where(cbf_valid, maps[f'{tissue}_cbf'], 0.0)
            for key in ('att', 'att_lower', 'att_upper'):
                maps[f'{tissue}_{key}'] = np.where(att_valid, maps[f'{tissue}_{key}'], 0.0)
            maps[f'{tissue}_kinstatus'] = status
            maps[f'{tissue}_valid'] = cbf_valid
            maps[f'{tissue}_att_valid'] = att_valid

        # The tissue curves are no longer needed; free them before the CSF regression.
        del fit['beta_gm'], fit['beta_wm']

        # CSF quality control: the signed least-squares CBF of the CSF curve against a
        # unit-CBF tissue curve (no arterial term), calibrated by the voxel's M0. CSF has no
        # tissue T1 in the tissue-T1 model, so the curve always decays with blood T1. Its ATT is the
        # voxel's uncorrected ATT where that was estimated, and otherwise the median of those,
        # so the value is defined in CSF voxels, where the uncorrected CBF is usually zero.
        csf_fit = local_regression_csf(
            signal, pv_gm, pv_wm, pv_csf, mask=mask, valid=quant_valid, **lr_kwargs
        )
        csf_sel = (
            np.isin(csf_fit['status'], (1, 2, 3, 5))
            & (pv_csf >= self.inputs.min_pv)
            & calibratable
        )
        unc_cbf = nb.load(self.inputs.mean_cbf).get_fdata()
        unc_att = nb.load(self.inputs.att).get_fdata()
        has_att = mask & np.isfinite(unc_cbf) & (unc_cbf > 0) & np.isfinite(unc_att)
        default_att = float(np.median(unc_att[has_att])) if has_att.any() else CSF_QC_DEFAULT_ATT
        basis_att = np.where(has_att, unc_att, default_att)
        n_sel = int(csf_sel.sum())
        basis = gkm_predict(
            np.ones(n_sel),
            basis_att[csf_sel],
            np.zeros(n_sel),
            np.zeros(n_sel),
            m0[csf_sel],
            plds[csf_sel],
            **kinetics,
        )
        csf_cbf_values, csf_kinres_values, defined = csf_kinetic_projection(
            csf_fit['beta_csf'][csf_sel], basis
        )
        csf_status = csf_fit['status'].copy()
        undefined = np.zeros(shape, bool)
        undefined[csf_sel] = ~defined
        csf_status[undefined] = 4
        csf_cbf = np.zeros(shape)
        csf_cbf[csf_sel] = csf_cbf_values
        csf_kinres = np.zeros(shape)
        csf_kinres[csf_sel] = csf_kinres_values
        csf_valid = csf_sel & ~undefined
        csf_summary = summarize_csf_qc(csf_cbf, csf_valid, pv_csf)

        outputs = {
            'gm_cbf': (maps['gm_cbf'], np.float32),
            'wm_cbf': (maps['wm_cbf'], np.float32),
            'gm_valid': (maps['gm_valid'], np.uint8),
            'wm_valid': (maps['wm_valid'], np.uint8),
            'gm_att': (maps['gm_att'], np.float32),
            'wm_att': (maps['wm_att'], np.float32),
            'gm_att_valid': (maps['gm_att_valid'], np.uint8),
            'wm_att_valid': (maps['wm_att_valid'], np.uint8),
            'gm_att_lower': (maps['gm_att_lower'], np.float32),
            'gm_att_upper': (maps['gm_att_upper'], np.float32),
            'wm_att_lower': (maps['wm_att_lower'], np.float32),
            'wm_att_upper': (maps['wm_att_upper'], np.float32),
            'gm_kinstatus': (maps['gm_kinstatus'], np.uint8),
            'wm_kinstatus': (maps['wm_kinstatus'], np.uint8),
            'gm_kinres': (maps['gm_kinres'], np.float32),
            'wm_kinres': (maps['wm_kinres'], np.float32),
            'status': (fit['status'], np.uint8),
            'condition': (fit['condition'], np.float32),
            'gm_noise': (fit['noise_gm'], np.float32),
            'wm_noise': (fit['noise_wm'], np.float32),
            'rmse': (fit['rmse'], np.float32),
            'csf_cbf': (csf_cbf, np.float32),
            'csf_status': (csf_status, np.uint8),
            'csf_noise': (np.where(csf_valid, csf_fit['noise_csf'], 0.0), np.float32),
            'csf_kinres': (csf_kinres, np.float32),
        }
        header = nb.load(self.inputs.m0_calib).header
        for name, (data, dtype) in outputs.items():
            out_img = nb.Nifti1Image(data.astype(dtype), affine, header)
            out_img.set_data_dtype(dtype)
            self._results[name] = fname_presuffix(
                self.inputs.deltam,
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
            'KineticModel': (
                'Each tissue curve (its delta-M amplitude at every delay, from local '
                'regression of each volume) was fit with the general kinetic model of the '
                'uncorrected multi-delay fit, with a non-negative arterial nuisance term.'
            ),
            'MaximumCBF': LR_CBF_MAX,
            'MaximumArterialBloodVolume': LR_ABV_MAX,
            'CSFQualityControlDefaultATT': default_att,
            'ResidualUnits': 'delta-M',
            'Calibration': (
                "Each tissue curve is calibrated by the centre voxel's M0, the same M0 the "
                'uncorrected CBF was quantified with.'
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
