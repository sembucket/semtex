#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Field file filter by HMB, patterned after an example from TA. Feb 2017.
#
# This utility is meant to replace the bash script c2w, which sets the
# scalar field c into the place of z-velocity component w (in fact
# this just requires changing the header), zeros all the other fields.
#
# Also, set the solution time and step in the header to zero.
#
# The input file is assumed to contain uvcp, output uvwp, where c-->w
# and all the other fields u, v, and p are set to zero.
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
    
    parser = argparse.ArgumentParser(description="Rename scalar c to w,"
                                     " zero other fields")
    parser.add_argument('infile', metavar='file', type=str,
                   help='the input file')
    args = parser.parse_args()

    ff = fieldfile.Fieldfile(args.infile, "r")
    data = ff.read()

    data_dict = {} 
    for i, field in enumerate(ff.fields):
        data_dict[field] = data[i]

    # -- Expect uvcp.
    
    if ff.hdr.fields != 'uvcp':
        sys.exit ("c2w.py: expected fields uvcp, found something else")
    
    # -- Zero everything but c.
    
    data_dict['u'].fill(0.0)
    data_dict['v'].fill(0.0)
    data_dict['p'].fill(0.0)
    
    # -- Write output.

    ff.hdr.fields = 'uvwp'
    ff.hdr.time   = 0.0
    ff.hdr.step   = 0
    ff_out = fieldfile.Fieldfile('/dev/stdout', 'w', ff.hdr)
    ff_out.write(data)
