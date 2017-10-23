#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# field file filter
# by Thomas Albrecht
#
# develop at /home/albrecht/tud/Semtex/semtex-work/utils

import numpy as np
import sys
from pdb import pm
import argparse
import fieldfile
import string
import copy
import logging

if __name__ == "__main__":
    """
    """
    
    parser = argparse.ArgumentParser(description="ff.py -- field file filter")
    parser.add_argument("-r", "--rename", help="rename fields, -r uAvb renames u->A and v->B")
    parser.add_argument("-K", "--keep", help="keep fields. This is executed after a possible --rename.")
    parser.add_argument("-s", "--scale", dest="scale", type=float, metavar='FACTOR', help="scale all fields by FACTOR")
    parser.add_argument("-n", "--normalise", action="store_true", help="normalise fields individually using maximum value")
    parser.add_argument("-m", "--minalise", action="store_true", help="normalise fields individually using minimum value")
    parser.add_argument("-v", "--velmag", action="store_true", help="add velocity magnitude")
    parser.add_argument('infile', metavar='file', type=str, nargs='+',
                   help='the file')
    args = parser.parse_args()

    ff = fieldfile.Fieldfile(args.infile[0], "r")
    data = ff.read()
    input_fields = ff.hdr.fields
    #elmt_wise = data.reshape((ff.nflds, ff.nz, ff.nel, ff.ns, ff.nr)) # works!

    if args.scale:
        data *= args.scale

    if args.normalise:
        for i, field in enumerate(ff.fields):
            data[i] /= max(data[i])

    if args.minalise:
        for i, field in enumerate(ff.fields):
            data[i] /= min(data[i])
            
    data_dict = {} 
    for i, field in enumerate(ff.fields):
        data_dict[field] = data[i]

    if args.velmag:
        data_dict['m'] = (data_dict['u']**2 + data_dict['v']**2 + data_dict['w']**2)**0.5
        input_fields += 'm'

    if not args.rename:
        renamed_fields = input_fields
    else:
        if len(args.rename) % 2 != 0:
            logging.error("rename arg must be pairs!")
            sys.exit(-1)
        src = args.rename[::2]
        dst = args.rename[1::2]
        for old_key in src:
            if old_key not in input_fields:
                logging.warn('field %s nominated for renaming not present.' % old_key)
        src2dst = dict(zip(src, dst))
        renamed_fields = ""
        for old_key in input_fields:
            if old_key in src2dst.keys():
                new_key = src2dst[old_key]
                renamed_fields += new_key
                data_dict[new_key] = data_dict.pop(old_key)
            else:
                renamed_fields += old_key

    if args.keep:
        renamed_fields = args.keep

       
    # -- write output
    hdr = copy.copy(ff.hdr)
    hdr.fields = ''.join(renamed_fields)
    ff_out = fieldfile.Fieldfile('/dev/stdout', 'w', hdr)
    ff_out.write(np.vstack(data_dict[f] for f in renamed_fields))
    
    
    #    print data
    #    - of nr * ns nodes
    #    - of 11th element
    #    - of first z-plane
    #    - of second field (in this case, v)
#    print elmt_wise[:,:,10,0,1]
