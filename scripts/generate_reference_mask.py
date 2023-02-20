#!/usr/bin/env python
import sys
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from aslprep.workflows.asl.util import init_asl_reference_wf


def sink_mask_file(in_file, orig_file, out_dir):
    import os
    from nipype.utils.filemanip import fname_presuffix, copyfile
    os.makedirs(out_dir, exist_ok=True)
    out_file = fname_presuffix(orig_file, suffix='_mask', newpath=out_dir)
    copyfile(in_file, out_file, copy=True, use_hardlink=True)
    return out_file


def init_main_wf(asl_file, out_dir, base_dir=None, name='main_wf'):
    wf = init_asl_reference_wf(omp_nthreads=4,
                                name=name)
    wf.base_dir = base_dir
    wf.inputs.inputnode.asl_file = asl_file

    sink = pe.Node(niu.Function(function=sink_mask_file),
                   name='sink')
    sink.inputs.out_dir = out_dir
    sink.inputs.orig_file = asl_file
    wf.connect([
        (wf.get_node('outputnode'), sink, [('asl_mask', 'in_file')]),
    ])
    return wf


def main():
    main_wf = init_main_wf(sys.argv[1], sys.argv[2])
    main_wf.run(plugin='MultiProc')


if __name__ == '__main__':
    main()
