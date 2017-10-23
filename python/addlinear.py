#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# Field file filter by HMB, patterned after an example from TA. Feb 2017.
#
# Expect to find an input file with fields uvwpf.  Set w += f and output uvwp.
# ----------------------------------------------------------------------------

import numpy as np
import sys
import argparse
import fieldfile
import string
import copy

if __name__ == "__main__":
    """
    """
    
    parser = argparse.ArgumentParser(description="Rename add trailing field f,"
                                     " to w")
    parser.add_argument('infile', metavar='file', type=str,
                   help='the input file')
    args = parser.parse_args()

    ff = fieldfile.Fieldfile(args.infile, "r")
    data = ff.read()

    data_dict = {} 
    for i, field in enumerate(ff.fields):
        data_dict[field] = data[i]

    # -- Expect uvwpf.
    
    if ff.hdr.fields != 'uvwpf':
        sys.exit ("addlinear.py: expected fields uvwpf, found something else")
    
    # -- Add f to w.
    
    data_dict['w'] += data_dict['f']
    
    # -- Write output.

    ff.hdr.fields = 'uvwp'
    ff.hdr.time   = 0.0
    ff.hdr.step   = 0
    ff_out = fieldfile.Fieldfile('/dev/stdout', 'w', ff.hdr)
    ff_out.write(data[0:4])
