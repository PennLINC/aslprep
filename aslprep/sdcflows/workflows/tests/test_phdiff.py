"""Test phase-difference type of fieldmaps."""
import pytest
from niworkflows.interfaces.bids import DerivativesDataSink
from nipype.pipeline import engine as pe

from ..phdiff import init_phdiff_wf, Workflow


@pytest.mark.parametrize('dataset', [
    'ds001600',
    'testdata',
])
def test_phdiff(bids_layouts, tmpdir, output_path, dataset, workdir):
    """Test creation of the workflow."""
    tmpdir.chdir()

    extra_entities = {}
    if dataset == 'ds001600':
        extra_entities['acquisition'] = 'v4'

    data = bids_layouts[dataset]
    wf = Workflow(name='phdiff_%s' % dataset)
    phdiff_wf = init_phdiff_wf(omp_nthreads=2)
    phdiff_wf.inputs.inputnode.magnitude = data.get(
        suffix=['magnitude1', 'magnitude2'],
        return_type='file',
        extension=['.nii', '.nii.gz'],
        **extra_entities)

    phdiff_files = data.get(
        suffix='phasediff',
        extension=['.nii', '.nii.gz'],
        **extra_entities)

    phdiff_wf.inputs.inputnode.phasediff = [
        (ph.path, ph.get_metadata()) for ph in phdiff_files]

    if output_path:
        from ...interfaces.reportlets import FieldmapReportlet
        rep = pe.Node(FieldmapReportlet(reference_label='Magnitude'), 'simple_report')
        rep.interface._always_run = True

        ds_report = pe.Node(DerivativesDataSink(
            base_directory=str(output_path), out_path_base='sdcflows', datatype='figures',
            suffix='fieldmap', desc='phasediff', dismiss_entities='fmap'), name='ds_report')
        ds_report.inputs.source_file = phdiff_files[0].path

        dsink_fmap = pe.Node(DerivativesDataSink(
            base_directory=str(output_path), dismiss_entities='fmap',
            desc='phasediff', suffix='fieldmap'), name='dsink_fmap')
        dsink_fmap.interface.out_path_base = 'sdcflows'
        dsink_fmap.inputs.source_file = phdiff_files[0].path

        wf.connect([
            (phdiff_wf, rep, [
                ('outputnode.fmap', 'fieldmap'),
                ('outputnode.fmap_ref', 'reference'),
                ('outputnode.fmap_mask', 'mask')]),
            (rep, ds_report, [('out_report', 'in_file')]),
            (phdiff_wf, dsink_fmap, [('outputnode.fmap', 'in_file')]),
        ])
    else:
        wf.add_nodes([phdiff_wf])

    if workdir:
        wf.base_dir = str(workdir)

    wf.run()


def test_phases(bids_layouts, tmpdir, output_path, workdir):
    """Test creation of the workflow."""
    tmpdir.chdir()

    data = bids_layouts['ds001600']
    wf = Workflow(name='phases_ds001600')
    phdiff_wf = init_phdiff_wf(omp_nthreads=2)
    phdiff_wf.inputs.inputnode.magnitude = data.get(
        suffix=['magnitude1', 'magnitude2'],
        acquisition='v2',
        return_type='file',
        extension=['.nii', '.nii.gz'])

    phdiff_files = data.get(suffix=['phase1', 'phase2'], acquisition='v2',
                            extension=['.nii', '.nii.gz'])

    phdiff_wf.inputs.inputnode.phasediff = [
        (ph.path, ph.get_metadata()) for ph in phdiff_files]

    if output_path:
        from ...interfaces.reportlets import FieldmapReportlet
        rep = pe.Node(FieldmapReportlet(reference_label='Magnitude'), 'simple_report')
        rep.interface._always_run = True

        ds_report = pe.Node(DerivativesDataSink(
            base_directory=str(output_path), out_path_base='sdcflows', datatype='figures',
            suffix='fieldmap', desc='twophases', dismiss_entities='fmap'), name='ds_report')
        ds_report.inputs.source_file = phdiff_files[0].path

        dsink_fmap = pe.Node(DerivativesDataSink(
            base_directory=str(output_path), suffix='fieldmap', desc='twophases',
            dismiss_entities='fmap'), name='dsink_fmap')
        dsink_fmap.interface.out_path_base = 'sdcflows'
        dsink_fmap.inputs.source_file = phdiff_files[0].path

        wf.connect([
            (phdiff_wf, rep, [
                ('outputnode.fmap', 'fieldmap'),
                ('outputnode.fmap_ref', 'reference'),
                ('outputnode.fmap_mask', 'mask')]),
            (rep, ds_report, [('out_report', 'in_file')]),
            (phdiff_wf, dsink_fmap, [('outputnode.fmap', 'in_file')]),
        ])
    else:
        wf.add_nodes([phdiff_wf])

    if workdir:
        wf.base_dir = str(workdir)

    wf.run()
