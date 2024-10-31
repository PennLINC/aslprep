# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Aggregate QC measures across all subjects in dataset."""
import os
from pathlib import Path

import pandas as pd


def get_parser():
    """Build parser object."""
    from argparse import ArgumentParser, RawTextHelpFormatter

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument('aslprep_dir', action='store', type=Path, help='aslprep output dir')

    parser.add_argument('output_prefix', action='store', type=str, help='output prefix for group')
    return parser


def main():
    """Run the workflow."""
    opts = get_parser().parse_args()

    allsubj_dir = os.path.abspath(opts.aslprep_dir)
    outputfile = os.getcwd() + '/' + str(opts.output_prefix) + '_allsubjects_qc.tsv'

    qclist = []
    for r, d, f in os.walk(allsubj_dir):
        for filex in f:
            if filex.endswith('desc-qualitycontrol_cbf.tsv'):
                qclist.append(r + '/' + filex)

    datax = pd.read_table(qclist[0])
    for i in range(1, len(qclist)):
        dy = pd.read_table(qclist[i])
        datax = pd.concat([datax, dy])

    datax.to_csv(outputfile, index=None, sep='\t')


if __name__ == '__main__':
    raise RuntimeError(
        'this should be run after running aslprep;\nit required installation of  aslprep'
    )
