"""Test pepolar type of fieldmaps."""
from os import cpu_count
import pytest
from niworkflows.interfaces.bids import DerivativesDataSink
from nipype.pipeline import engine as pe

from ..pepolar import (
    check_pes, _split_epi_lists, init_prepare_epi_wf, init_pepolar_unwarp_wf
)


def test_split_epi_lists(bids_layouts, tmpdir):
    """Test preparation workflow."""
    tmpdir.chdir()

    layout = bids_layouts['testdata']

    bold = layout.get(suffix='bold', direction='LR',
                      extension=['.nii.gz', '.nii'])[0]
    epidata = layout.get(suffix='epi', desc=None, extension=['.nii.gz', '.nii'])

    # EPI fmaps in HCP101006 have 3 volumes each
    a, b = _split_epi_lists(
        in_files=[(im.path, im.get_metadata()['PhaseEncodingDirection'])
                  for im in epidata],
        pe_dir=bold.get_metadata()['PhaseEncodingDirection']
    )

    assert len(a) == 3
    assert len(b) == 3

    # If we append the BOLD aligned (not that you would do this), then the
    # second array should have 53 volumes (50 from _bold, 3 from _epi)
    a, b = _split_epi_lists(
        in_files=[(im.path, im.get_metadata()['PhaseEncodingDirection'])
                  for im in epidata + [bold]],
        pe_dir=bold.get_metadata()['PhaseEncodingDirection']
    )

    assert len(a) == 3
    assert len(b) == 53


def test_prepare_epi_wf0(bids_layouts, tmpdir):
    """Test preparation workflow."""
    tmpdir.chdir()

    layout = bids_layouts['testdata']

    bold = layout.get(suffix='bold', direction='LR',
                      extension=['.nii.gz', '.nii'])[0]
    epidata = layout.get(suffix='epi', direction='LR',
                         desc=None, extension=['.nii.gz', '.nii'])

    with pytest.raises(ValueError):
        check_pes([
            (im.path, im.get_metadata()['PhaseEncodingDirection'])
            for im in epidata], bold.get_metadata()['PhaseEncodingDirection'])


def test_prepare_epi_wf1(bids_layouts, tmpdir):
    """Test preparation workflow."""
    tmpdir.chdir()

    layout = bids_layouts['testdata']

    bold = layout.get(suffix='bold', direction='LR',
                      extension=['.nii.gz', '.nii'])[0]
    boldref = layout.get(suffix='boldref', direction='LR',
                         extension=['.nii.gz', '.nii'])[0]
    epidata = layout.get(suffix='epi', desc=None, extension=['.nii.gz', '.nii'])

    matched_pe = check_pes([(im.path, im.get_metadata()['PhaseEncodingDirection'])
                            for im in epidata], bold.get_metadata()['PhaseEncodingDirection'])

    assert matched_pe is True

    wf = init_prepare_epi_wf(omp_nthreads=1, matched_pe=matched_pe)
    wf.inputs.inputnode.maps_pe = [(im.path, im.get_metadata()['PhaseEncodingDirection'])
                                   for im in epidata]
    wf.inputs.inputnode.ref_brain = boldref.path
    wf.inputs.inputnode.epi_pe = bold.get_metadata()['PhaseEncodingDirection']


def test_prepare_epi_wf2(bids_layouts, tmpdir):
    """Test preparation workflow."""
    tmpdir.chdir()

    layout = bids_layouts['testdata']

    bold = layout.get(suffix='bold', direction='LR',
                      extension=['.nii.gz', '.nii'])[0]
    boldref = layout.get(suffix='boldref', direction='LR',
                         extension=['.nii.gz', '.nii'])[0]
    epidata = layout.get(suffix='epi', direction='RL',
                         desc=None, extension=['.nii.gz', '.nii'])

    matched_pe = check_pes([(im.path, im.get_metadata()['PhaseEncodingDirection'])
                            for im in epidata], bold.get_metadata()['PhaseEncodingDirection'])

    assert matched_pe is False

    wf = init_prepare_epi_wf(omp_nthreads=1, matched_pe=matched_pe)
    wf.inputs.inputnode.maps_pe = [(im.path, im.get_metadata()['PhaseEncodingDirection'])
                                   for im in epidata]
    wf.inputs.inputnode.ref_brain = boldref.path
    wf.inputs.inputnode.epi_pe = bold.get_metadata()['PhaseEncodingDirection']


@pytest.mark.parametrize('dataset', [
    'ds001600',
    'testdata',
])
def test_pepolar_wf1(bids_layouts, output_path, dataset, workdir):
    """Test preparation workflow."""
    layout = bids_layouts[dataset]

    if dataset == 'testdata':
        bold = layout.get(suffix='bold', direction='LR',
                          extension=['.nii.gz', '.nii'])[0]
        boldref = layout.get(suffix='boldref', direction='LR', desc='brain',
                             extension=['.nii.gz', '.nii'])[0]
        epidata = layout.get(suffix='epi', desc=None, extension=['.nii.gz', '.nii'])
    elif dataset == 'ds001600':
        bold = layout.get(suffix='bold', acquisition='AP',
                          extension=['.nii.gz', '.nii'])[0]
        epidata = layout.get(suffix='epi', extension=['.nii.gz', '.nii'])

    matched_pe = check_pes([(im.path, im.get_metadata()['PhaseEncodingDirection'])
                            for im in epidata], bold.get_metadata()['PhaseEncodingDirection'])

    wf = init_pepolar_unwarp_wf(omp_nthreads=cpu_count(), matched_pe=matched_pe)
    wf.inputs.inputnode.fmaps_epi = [(im.path, im.get_metadata()['PhaseEncodingDirection'])
                                     for im in epidata]
    wf.inputs.inputnode.epi_pe_dir = bold.get_metadata()['PhaseEncodingDirection']

    if output_path:
        from nipype.interfaces import utility as niu
        from ..pepolar import Workflow
        from ...interfaces.reportlets import FieldmapReportlet

        boiler = Workflow(name='pepolar_%s' % dataset)

        split_field = pe.Node(niu.Function(function=_split_field), name='split_field')

        if dataset == 'ds001600':
            from niworkflows.func.util import init_bold_reference_wf
            gen_ref = init_bold_reference_wf(
                omp_nthreads=cpu_count(),
                bold_file=bold.path)
            boiler.connect([
                (gen_ref, wf, [
                    ('outputnode.ref_image', 'inputnode.in_reference'),
                    ('outputnode.ref_image_brain', 'inputnode.in_reference_brain')])
            ])
        else:
            wf.inputs.inputnode.in_reference_brain = boldref.path
            wf.inputs.inputnode.in_reference = boldref.path

        rep = pe.Node(FieldmapReportlet(reference_label='EPI Reference'), 'simple_report')
        rep.interface._always_run = True
        dsink = pe.Node(DerivativesDataSink(
            base_directory=str(output_path), out_path_base='sdcflows', datatype="figures",
            suffix='fieldmap', desc='pepolar', dismiss_entities='fmap'), name='dsink')
        dsink.inputs.source_file = epidata[0].path

        boiler.connect([
            (wf, split_field, [
                ('inputnode.epi_pe_dir', 'pe_dir'),
                ('outputnode.out_warp', 'in_field')]),
            (split_field, rep, [
                ('out', 'fieldmap')]),
            (wf, rep, [
                # ('outputnode.out_warp', 'fieldmap'),
                ('outputnode.out_reference_brain', 'reference'),
                ('outputnode.out_mask', 'mask')]),
            (rep, dsink, [('out_report', 'in_file')]),
        ])

        if workdir:
            boiler.base_dir = str(workdir)
        boiler.run(plugin='MultiProc', plugin_args={'n_proc': cpu_count()})


def _split_field(in_field, pe_dir):
    from os.path import abspath
    import numpy as np
    import nibabel as nb
    axis = 'ijk'.index(pe_dir[0])
    im = nb.load(in_field)
    data = np.squeeze(im.get_fdata())[..., axis]
    dirnii = nb.Nifti1Image(data, im.affine, im.header)
    dirnii.to_filename('fieldmap.nii.gz')
    return abspath('fieldmap.nii.gz')
