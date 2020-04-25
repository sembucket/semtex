#!/usr/bin/env python
#
# Interactive field file composition utility to combine two field
# files, with field selection.  Alternatively, given a single field
# file, choose (and perhaps rename) selected fields.  Write binary
# file (asked for name).
#
# 1. compose.py file.fld
# 2. compose.py file1.fld file2.fld
#
# Input files are assumed to be in a machine-compatible binary format.
# In the second case, the field files must conform (nr, ns, nel, nz).
#
# This utility cannot reorder output fields, and unless you rename/choose
# output field names appropriately, it possible to repeat field names.
#
# $Id: compose.py,v 9.3 2020/01/23 03:15:43 hmb Exp $
# ----------------------------------------------------------------------------

import numpy as np
import sys
import argparse
import fieldfile
import re

if (len(sys.argv) == 2):
    mode = 1
    file1 = fieldfile.Fieldfile(sys.argv[1], "r")
elif (len(sys.argv) == 3):
    mode = 2
    file1 = fieldfile.Fieldfile(sys.argv[1], "r")
    file2 = fieldfile.Fieldfile(sys.argv[2], "r")
    if ((file1.hdr.geometry.nr  != file2.hdr.geometry.nr)   or
        (file1.hdr.geometry.ns  != file2.hdr.geometry.ns)   or
        (file1.hdr.geometry.nz  != file2.hdr.geometry.nz)   or
        (file1.hdr.geometry.nel != file2.hdr.geometry.nel)) :
        print "The two input files do not conform"
        sys.exit(1)
else:
    print "Usage: compose.py file.fld [anotherfile.fld]"
    sys.exit(1)

if (mode == 1):
    print 'Input has these fields, supply a string below to',
    print 'indicate/rename ones you want:'
    print file1.hdr.fields
    required_fields = raw_input()
    if (len(required_fields) < len(file1.hdr.fields)):
        print required_fields, ': cannot shorter than :', file1.hdr.fields
        sys.exit(1)
    wanted = ''.join(re.findall(r'\S+', required_fields))
    
    outfile = raw_input("type in an output file name: ")
    newhdr = file1.hdr
    newhdr.fields = wanted
    ofile = fieldfile.Fieldfile(outfile, "w", newhdr)

    ofile.data = np.zeros((len(ofile.hdr.fields),ofile.ntot))

    j = 0
    for i, field in enumerate(file1.fields):
        buf = file1.f.read(file1.ntot*8)
        if (required_fields[i].isspace()) == False:
            ofile.data[j] = np.fromstring(buf, np.float64, -1)
            j += 1

    ofile.write(ofile.data)
    file1.close()
    ofile.close()

else:
    print 'Input1 has these fields, supply a string below to',
    print 'indicate/rename ones you want:'
    print file1.hdr.fields
    required_fields1 = raw_input()
    if (len(required_fields1) < len(file1.hdr.fields)):
        print required_fields1, ': cannot be shorter than :', file1.hdr.fields
        sys.exit(1)
    wanted1 = ''.join(re.findall(r'\S+', required_fields1))
    
    print 'Input2 has these fields, supply a string below to',
    print 'indicate/rename ones you want:'
    print file2.hdr.fields
    required_fields2 = raw_input()
    if (len(required_fields2) < len(file2.hdr.fields)):
        print required_fields2, ': cannot be shorter than :', file2.hdr.fields
        sys.exit(1)    
    wanted2 = ''.join(re.findall(r'\S+', required_fields2))

    outfile = raw_input("type in an output file name: ")
    newhdr  = file1.hdr
    newhdr.fields = wanted1 + wanted2
    print 'new file will have fields: ', newhdr.fields
    ofile = fieldfile.Fieldfile(outfile, 'w', newhdr)

    ofile.data = np.zeros((len(ofile.hdr.fields),ofile.ntot))

    j = 0
    for i, field in enumerate(file1.fields):
        buf = file1.f.read(file1.ntot*8)
        if (required_fields1[i].isspace()) == False:
            ofile.data[j] = np.fromstring(buf, np.float64, -1)
            j += 1
    for i, field in enumerate(file2.fields):
        buf = file2.f.read(file2.ntot*8)
        if (required_fields2[i].isspace()) == False:
            ofile.data[j] = np.fromstring(buf, np.float64, -1)
            j += 1

    ofile.write(ofile.data)            
    file1.close()
    file2.close()
    ofile.close()
