#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
"""
zero_plane -- zero all but given planes/modes
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
import logging

def range_to_list(rng):
    """accept a string "x-y" or "x". return a list x, x+1, ... y"""
    if '-' in rng:
        x0, x1 = [int(item) for item in rng.split('-')]
        return range(x0, x1)
    else:
        return [int(rng)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="zero_plane -- zero all but given planes/modes")
    parser.add_argument("-z", help="list of planes, e.g. 2 or 0,1,7-9. A range excludes the upper boundary: 7-9 -> 7,8", default=None)
    parser.add_argument("-m", help="list of modes", default=None)
    parser.add_argument("-x", "--extract", action="store_true", help="Extract the given planes/modes (similar to xplane).")
    parser.add_argument('infile', metavar='file', type=str, nargs='+',
                   help='the file')
    args = parser.parse_args()
            
    
    Z = []
    if args.z:
        for item in args.z.split(','):
            Z += range_to_list(item)
            
    M = []     
    if args.m:
        for item in args.m.split(','):
            M += range_to_list(item)
        
        Z += [2*the_mode   for the_mode in M]
        Z += [2*the_mode+1 for the_mode in M]
        Z.sort()

    if len(Z) < 1:
        logging.error("No planes given.")
        sys.exit(-1)
        

    ff = fieldfile.Fieldfile(args.infile[0], "r")
    data = ff.read()
    elmt_wise = data.reshape((ff.nflds, ff.nz, ff.nel, ff.ns, ff.nr))

    if args.extract:
        z_wise = elmt_wise.reshape(ff.nflds, ff.nz, ff.nxy)
        mask = z_wise < -1e99 # create an all-false mask
        for z in range(ff.nz):
            if z in Z:
                mask[:,z,:] = True
        elmt_wise = z_wise[mask]
        ff.hdr.geometry.nz = len(Z)
    else:
        for z in range(ff.nz):
            if z not in Z:
                elmt_wise[:, z, :, :, :] = 0.

    ff_out = fieldfile.Fieldfile('/dev/stdout', 'w', ff.hdr)
    ff_out.write(elmt_wise.flatten())
