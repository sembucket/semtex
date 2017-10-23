#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
print min max per field per plane
"""
#
# by Thomas Albrecht
#

import numpy as np
import sys
from pdb import pm
import argparse
import fieldfile
import string
import copy

def fourier_name(z):
    if (z % 2) == 0:
        return "%3i Re" % (z/2)
    else:
        return "%3i Im" % (z/2)

if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser(description="fpminmax.py -- print min max per field per plane")
    parser.add_argument("-e", "--eps", type=float, help="print 0 if value below EPS", default=1e-15)
    parser.add_argument('infile', metavar='file', type=str, nargs='+',
                   help='the file')
    args = parser.parse_args()

    ff = fieldfile.Fieldfile(args.infile[0], "r")
    data = ff.read()
    elmt_wise = data.reshape((ff.nflds, ff.nz, ff.nel, ff.ns, ff.nr))
    
    if np.isnan(data).any():
        print "NaN"
    if np.isinf(data).any():
        print "Inf"
    
    print "#  z     FM ",
    for field in ff.fields:
        print "%9s_min %8s_max " % (field, field),
    print
    for z in range(ff.nz):
        print "%4i %s" % (z, fourier_name(z)),
        for field in range(ff.nflds):
            fmin = elmt_wise[field,z].flatten().min()
            fmax = elmt_wise[field,z].flatten().max()
            if abs(fmin) < args.eps: fmin = 0.
            if abs(fmax) < args.eps: fmax = 0.
            print "  % 12g % 12g" % (fmin, fmax),
        print

