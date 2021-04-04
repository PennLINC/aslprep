# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""aggregate qc of all the subjects"""
import os
import glob as glob 
import pandas as pd 

def get_parser():
    """Build parser object"""
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument('aslprep_dir', action='store', type=Path, help='aslprep output dir')
    
    parser.add_argument('output_prefix', action='store', type=Path, help='output prefix for group')
    return parser


def main():
    """Entry point"""
    
    opts = get_parser().parse_args()

    allsubj_dir = opts.aslprep_dir
    outputfile  = opts.output_prefix + '_allsubjects_qc.csv'
    
    qclist=[]
    for r, d, f in os.walk(allsubj_dir ):
        for filex in f:
            if filex.endswith("_quality_control_cbf.csv"):
                qclist.append(filex)
    datax = []
    for i in qclist:
        dy = pd.read_csv(i)
        datax.append(dy)
    
    df = pd.concat(datax)   
    
    df.to_csv(outputfile,index=None)

if __name__ == '__main__':
    raise RuntimeError("aslprep/cli/run.py should not be run directly;\n"
                       "Please `pip install` aslprep and use the `aslprep` command")
