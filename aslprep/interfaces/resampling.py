# Adapted from fmriprep.interfaces.resampling (fMRIPrep 25.2),
# Copyright The NiPreps Developers, licensed under the Apache License, Version 2.0.
# See the "External code" section of LICENSE.md.
"""Single-shot resampling with gradient nonlinearity correction."""

import nibabel as nb
import nitransforms as nt
import numpy as np
from fmriprep.interfaces.resampling import ResampleSeries as _FMRIPrepResampleSeries
from fmriprep.interfaces.resampling import ResampleSeriesInputSpec as _FMRIPrepInputSpec
from fmriprep.interfaces.resampling import resample_series, resample_vol
from fmriprep.utils.transforms import load_transforms
from nipype.interfaces.base import File, isdefined, traits
from nipype.utils.filemanip import fname_presuffix
from scipy import ndimage as ndi
from sdcflows.utils.tools import ensure_positive_cosines

#: Converts ITK's LPS displacement vectors to RAS, and back.
_LPS_TO_RAS = np.array([-1.0, -1.0, 1.0])


class _ResampleSeriesInputSpec(_FMRIPrepInputSpec):
    gradwarp_field = File(
        exists=True,
        desc=(
            'ITK displacement field (LPS, mm) for gradient nonlinearity correction, '
            "as written by TORTOISE's CreateNonlinearityDisplacementMap. "
            'It maps a point in the corrected image to where it was recorded.'
        ),
    )
    gradwarp_jacobian = traits.Bool(
        True,
        usedefault=True,
        desc='Whether to modulate intensities by the Jacobian determinant of gradwarp_field',
    )


class ResampleSeries(_FMRIPrepResampleSeries):
    """Resample a time series, applying gradient nonlinearity, susceptibility distortion,
    and motion correction simultaneously.

    Without ``gradwarp_field``, this is fMRIPrep's ``ResampleSeries``.
    """

    input_spec = _ResampleSeriesInputSpec

    def _run_interface(self, runtime):
        """Resample ``in_file``, including the gradient field if one is given.

        Parameters
        ----------
        runtime : nipype.interfaces.base.support.Bunch
            Nipype runtime object.

        Returns
        -------
        runtime : nipype.interfaces.base.support.Bunch
            The same runtime object. ``out_file`` is set in the results.
        """
        if not isdefined(self.inputs.gradwarp_field):
            return super()._run_interface(runtime)

        out_path = fname_presuffix(self.inputs.in_file, suffix='resampled', newpath=runtime.cwd)

        source = nb.load(self.inputs.in_file)
        target = nb.load(self.inputs.ref_file)
        fieldmap = nb.load(self.inputs.fieldmap) if self.inputs.fieldmap else None

        nvols = source.shape[3] if source.ndim > 3 else 1

        # No transforms appear Undefined, pass as empty list
        transforms = load_transforms(self.inputs.transforms or [], self.inputs.inverse)

        pe_dir = self.inputs.pe_dir
        ro_time = self.inputs.ro_time
        pe_info = None

        if pe_dir and ro_time:
            pe_axis = 'ijk'.index(pe_dir[0])
            pe_flip = pe_dir.endswith('-')

            # Nitransforms displacements are positive
            source, axcodes = ensure_positive_cosines(source)
            axis_flip = axcodes[pe_axis] in 'LPI'

            pe_info = [(pe_axis, -ro_time if (axis_flip ^ pe_flip) else ro_time)] * nvols

        resampled = resample_image(
            source=source,
            target=target,
            transforms=transforms,
            fieldmap=fieldmap,
            pe_info=pe_info,
            gradwarp=GradwarpField(nb.load(self.inputs.gradwarp_field)),
            jacobian=self.inputs.jacobian,
            gradwarp_jacobian=self.inputs.gradwarp_jacobian,
            nthreads=self.inputs.num_threads,
            output_dtype=self.inputs.output_data_type,
            order=self.inputs.order,
            mode=self.inputs.mode,
            cval=self.inputs.cval,
            prefilter=self.inputs.prefilter,
        )
        resampled.to_filename(out_path)

        self._results['out_file'] = out_path
        return runtime


class GradwarpField:
    """A gradient nonlinearity displacement field, evaluated in world (RAS) coordinates.

    Parameters
    ----------
    img : nibabel.Nifti1Image
        ITK displacement field: shape (X, Y, Z, 1, 3) or (X, Y, Z, 3), in LPS millimetres.

    Attributes
    ----------
    deltas : numpy.ndarray of shape (X, Y, Z, 3)
        The displacements, in RAS millimetres.
    affine : numpy.ndarray of shape (4, 4)
        The field grid's voxel-to-world (RAS) affine.

    Raises
    ------
    ValueError
        If ``img`` does not have the shape of an ITK displacement field.

    Notes
    -----
    The field is smooth (a low-order spherical harmonic expansion), so it is sampled linearly.
    Points beyond the field's grid take the displacement at the nearest edge of the grid.
    The field is defined on the ASL reference grid, so such points carry no ASL signal anyway,
    but this keeps target grids that extend past the ASL field of view (e.g., T1w or template
    space) free of the NaNs nitransforms' ``DenseFieldTransform`` would produce there.
    """

    def __init__(self, img: nb.Nifti1Image):
        data = np.asarray(img.dataobj, dtype='float64')
        if data.shape[-1] != 3 or data.ndim not in (4, 5):
            raise ValueError(
                f'Expected an ITK displacement field of shape (X, Y, Z, 1, 3), got {data.shape}.'
            )
        self.deltas = data.reshape(data.shape[:3] + (3,)) * _LPS_TO_RAS
        self.affine = img.affine
        self._ras2vox = np.linalg.inv(img.affine)
        self._jacobian = None

    def _index(self, points: np.ndarray) -> np.ndarray:
        """Convert world points to voxel indices on the field grid.

        Parameters
        ----------
        points : numpy.ndarray of shape (N, 3)
            Points in RAS world coordinates (mm).

        Returns
        -------
        ijk : numpy.ndarray of shape (3, N)
            Continuous voxel indices on the field grid.
        """
        return self._ras2vox[:3, :3] @ points.T + self._ras2vox[:3, 3:4]

    def map(self, points: np.ndarray) -> np.ndarray:
        """Map points in the corrected image to where they were recorded.

        Parameters
        ----------
        points : numpy.ndarray of shape (N, 3)
            Points in RAS world (scanner) coordinates (mm), in the corrected image.

        Returns
        -------
        recorded : numpy.ndarray of shape (N, 3)
            ``points`` plus the linearly interpolated displacement: where each point was
            recorded in the distorted (raw) image, in RAS world coordinates (mm).
        """
        ijk = self._index(points)
        displacements = np.stack(
            [
                ndi.map_coordinates(self.deltas[..., i], ijk, order=1, mode='nearest')
                for i in range(3)
            ],
            axis=1,
        )
        return points + displacements

    def jacobian(self) -> np.ndarray:
        r"""Compute the Jacobian determinant of :math:`x \mapsto x + u(x)` on the field grid.

        The result is computed once and cached.

        Returns
        -------
        det : numpy.ndarray of shape (X, Y, Z)
            The determinant of :math:`I + \partial u / \partial x` at each voxel of the field
            grid, with derivatives taken by central differences in world coordinates.
            Values above 1 mean that a corrected voxel draws from a larger raw volume.
        """
        if self._jacobian is not None:
            return self._jacobian
        # du_i/dj_k (derivatives along voxel axes), then the chain rule to world coordinates.
        grad_vox = np.stack(
            [np.stack(np.gradient(self.deltas[..., i]), axis=-1) for i in range(3)],
            axis=-2,
        )
        grad_world = grad_vox @ self._ras2vox[:3, :3]
        self._jacobian = np.linalg.det(np.eye(3) + grad_world)
        return self._jacobian

    def sample_jacobian(self, points: np.ndarray) -> np.ndarray:
        """Sample the Jacobian determinant at points in the corrected image.

        Parameters
        ----------
        points : numpy.ndarray of shape (N, 3)
            Points in RAS world (scanner) coordinates (mm), in the corrected image.

        Returns
        -------
        det : numpy.ndarray of shape (N,)
            The linearly interpolated Jacobian determinant (see :meth:`jacobian`),
            taking the nearest edge value beyond the field grid.
        """
        return ndi.map_coordinates(self.jacobian(), self._index(points), order=1, mode='nearest')


def resample_image(
    source: nb.Nifti1Image,
    target: nb.Nifti1Image,
    transforms: nt.TransformChain,
    fieldmap: nb.Nifti1Image | None,
    pe_info: list[tuple[int, float]] | None,
    gradwarp: GradwarpField | None = None,
    jacobian: bool = True,
    gradwarp_jacobian: bool = True,
    nthreads: int = 1,
    output_dtype: np.dtype | str | None = 'f4',
    order: int = 3,
    mode: str = 'constant',
    cval: float = 0.0,
    prefilter: bool = True,
) -> nb.Nifti1Image:
    """Resample a 3- or 4D image into a target space, applying gradient nonlinearity,
    head-motion, and susceptibility-distortion correction simultaneously.

    This is :func:`fmriprep.interfaces.resampling.resample_image` with one added step.
    Target coordinates are mapped through ``transforms`` (except head motion) into the world
    frame of the reference image, and through each volume's head-motion transform into that
    volume's world frame, which is where the scanner-fixed ``gradwarp`` field is evaluated.
    Only then are they mapped into source voxels, where the fieldmap shift is applied as in
    fMRIPrep. Unlike the susceptibility field, which moves with the head, the gradient field
    does not, so it must be evaluated after head motion, separately for each volume.

    Parameters
    ----------
    source : nibabel.Nifti1Image
        The 3D image or 4D series to resample.
    target : nibabel.Nifti1Image
        An image sampled in the target space.
    transforms : nitransforms.TransformChain
        Transforms mapping the target space to the source. If head-motion transforms
        (a :class:`nitransforms.linear.LinearTransformsMapping`) are included, they must
        be last, and map the reference image's world frame to each volume's.
    fieldmap : nibabel.Nifti1Image or None
        The fieldmap, in Hz, sampled in the target space.
    pe_info : list of tuple of (int, float), or None
        For each volume, the phase-encoding axis and signed readout time, as in fMRIPrep.
    gradwarp : GradwarpField or None, optional
        Gradient nonlinearity displacement field. If None, no gradient correction is applied
        and the result is that of fMRIPrep's ``resample_image``.
    jacobian : bool, optional
        Whether to modulate intensities by the Jacobian of the fieldmap shift.
    gradwarp_jacobian : bool, optional
        Whether to modulate intensities by the Jacobian determinant of ``gradwarp``.
    nthreads : int, optional
        Number of threads to use.
    output_dtype : numpy.dtype or str or None, optional
        Data type of the resampled array.
    order : int, optional
        Order of spline interpolation (default: 3).
    mode : str, optional
        How data are extended beyond their boundaries,
        as in :func:`scipy.ndimage.map_coordinates`.
    cval : float, optional
        Value used beyond the boundaries when ``mode`` is ``'constant'``.
    prefilter : bool, optional
        Whether to spline-prefilter the data when ``order`` > 1.

    Returns
    -------
    resampled_img : nibabel.Nifti1Image
        The source resampled into the target space.

    Raises
    ------
    ValueError
        If head-motion transforms are not last in ``transforms``.
    """
    if not isinstance(transforms, nt.TransformChain):
        transforms = nt.TransformChain([transforms])
    if isinstance(transforms[-1], nt.linear.LinearTransformsMapping):
        transform_list, hmc = list(transforms[:-1]), transforms[-1]
    else:
        if any(isinstance(xfm, nt.linear.LinearTransformsMapping) for xfm in transforms):
            classes = [xfm.__class__.__name__ for xfm in transforms]
            raise ValueError(f'HMC transforms must come last. Found sequence: {classes}')
        transform_list = list(transforms.transforms)
        hmc = []

    # Retrieve the RAS coordinates of the target space
    coordinates = nt.base.SpatialReference.factory(target).ndcoords.astype('f4')

    # We will operate in voxel space, so get the source affine
    vox2ras = source.affine
    ras2vox = np.linalg.inv(vox2ras)

    # Some identities to reduce special casing downstream
    if fieldmap is None:
        fieldmap = nb.Nifti1Image(np.zeros(target.shape[:3], dtype='f4'), target.affine)
    if pe_info is None:
        pe_info = [[0, 0] for _ in range(source.shape[-1])]

    if gradwarp is None:
        # fMRIPrep's path: one set of coordinates, with head motion applied in voxel space
        ref2vox = nt.TransformChain(transform_list + [nt.Affine(ras2vox)])
        resampled_data = resample_series(
            data=source.get_fdata(dtype='f4'),
            coordinates=ref2vox.map(coordinates).T.reshape((3, *target.shape[:3])),
            pe_info=pe_info,
            jacobian=jacobian,
            hmc_xfms=[ras2vox @ xfm.matrix @ vox2ras for xfm in hmc],
            fmap_hz=fieldmap.get_fdata(dtype='f4'),
            output_dtype=output_dtype,
            nthreads=nthreads,
            order=order,
            mode=mode,
            cval=cval,
            prefilter=prefilter,
        )
    else:
        reference_coordinates = (
            nt.TransformChain(transform_list).map(coordinates) if transform_list else coordinates
        )
        resampled_data = _resample_series_gradwarp(
            data=source.get_fdata(dtype='f4'),
            reference_coordinates=reference_coordinates,
            target_shape=target.shape[:3],
            ras2vox=ras2vox,
            hmc_matrices=[xfm.matrix for xfm in hmc],
            gradwarp=gradwarp,
            gradwarp_jacobian=gradwarp_jacobian,
            pe_info=pe_info,
            jacobian=jacobian,
            fmap_hz=fieldmap.get_fdata(dtype='f4'),
            output_dtype=output_dtype,
            nthreads=nthreads,
            order=order,
            mode=mode,
            cval=cval,
            prefilter=prefilter,
        )

    resampled_img = nb.Nifti1Image(resampled_data, target.affine, target.header)
    resampled_img.set_data_dtype('f4')
    # Preserve zooms of additional dimensions
    resampled_img.header.set_zooms(target.header.get_zooms()[:3] + source.header.get_zooms()[3:])

    return resampled_img


def _resample_series_gradwarp(
    *,
    data: np.ndarray,
    reference_coordinates: np.ndarray,
    target_shape: tuple[int, int, int],
    ras2vox: np.ndarray,
    hmc_matrices: list[np.ndarray],
    gradwarp: GradwarpField,
    gradwarp_jacobian: bool,
    pe_info: list[tuple[int, float]],
    jacobian: bool,
    fmap_hz: np.ndarray,
    output_dtype,
    nthreads: int,
    order: int,
    mode: str,
    cval: float,
    prefilter: bool,
) -> np.ndarray:
    """Resample each volume with the gradient field evaluated after its head motion.

    Parameters
    ----------
    data : numpy.ndarray of shape (X, Y, Z) or (X, Y, Z, T)
        The source data.
    reference_coordinates : numpy.ndarray of shape (N, 3)
        The target voxels' positions in the reference image's RAS world frame (mm),
        where N is the number of target voxels.
    target_shape : tuple of int
        Shape of the target grid, whose product is N.
    ras2vox : numpy.ndarray of shape (4, 4)
        World-to-voxel affine of the source data.
    hmc_matrices : list of numpy.ndarray of shape (4, 4)
        For each volume, the RAS-to-RAS affine from the reference image's world frame to
        that volume's. Empty for no head-motion correction.
    gradwarp : GradwarpField
        Gradient nonlinearity displacement field.
    gradwarp_jacobian : bool
        Whether to modulate intensities by the Jacobian determinant of ``gradwarp``,
        evaluated where each volume's tissue was in the scanner.
    pe_info : list of tuple of (int, float)
        For each volume, the phase-encoding axis and signed readout time.
    jacobian : bool
        Whether to modulate intensities by the Jacobian of the fieldmap shift.
    fmap_hz : numpy.ndarray of shape ``target_shape``
        The fieldmap, in Hz, sampled in the target space.
    output_dtype : numpy.dtype or str
        Data type of the resampled array.
    nthreads : int
        Number of volumes to resample in parallel.
    order : int
        Order of spline interpolation.
    mode : str
        How data are extended beyond their boundaries,
        as in :func:`scipy.ndimage.map_coordinates`.
    cval : float
        Value used beyond the boundaries when ``mode`` is ``'constant'``.
    prefilter : bool
        Whether to spline-prefilter the data when ``order`` > 1.

    Returns
    -------
    resampled : numpy.ndarray
        The resampled data, of shape ``target_shape`` for 3D ``data``,
        or ``target_shape + (T,)`` for 4D ``data``.
    """
    from concurrent.futures import ThreadPoolExecutor

    volumes = data[..., np.newaxis] if data.ndim == 3 else data
    nvols = volumes.shape[-1]
    out_array = np.zeros(target_shape + (nvols,), dtype=output_dtype, order='F')

    def _resample_one(volid):
        """Resample one volume into ``out_array``.

        Parameters
        ----------
        volid : int
            Index of the volume to resample.
        """
        points = reference_coordinates
        if hmc_matrices:
            points = nb.affines.apply_affine(hmc_matrices[volid], points)
        voxels = nb.affines.apply_affine(ras2vox, gradwarp.map(points))
        resample_vol(
            data=volumes[..., volid],
            coordinates=voxels.T.reshape((3, *target_shape)),
            pe_info=pe_info[volid],
            jacobian=jacobian,
            hmc_xfm=None,
            fmap_hz=fmap_hz,
            output=out_array[..., volid],
            order=order,
            mode=mode,
            cval=cval,
            prefilter=prefilter,
        )
        if gradwarp_jacobian:
            out_array[..., volid] *= gradwarp.sample_jacobian(points).reshape(target_shape)

    with ThreadPoolExecutor(max_workers=max(nthreads, 1)) as executor:
        list(executor.map(_resample_one, range(nvols)))

    return out_array[..., 0] if data.ndim == 3 else out_array
