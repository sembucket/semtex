#!/usr/bin/env python
#
# Usage: roll+streak.py mesh_file roll.dat streak.dat Rey
# 
# Combine u, v, p in roll file with w in streak file in order to
# generate a 2D/3C roll + streak flow.  The final argument, scale, is
# used to scale the roll variables relative to the steak variables;
# typically this is 1/Re.  The (meshpr-type) mesh file is used to
# generate structure information that gets used to make a semtex-style
# header.
# 
# $Id$
##############################################################################

import sys, os

try:
    mfilename = sys.argv[1]
    rfilename = sys.argv[2]
    sfilename = sys.argv[3]
    scale = 1.0/float(sys.argv[4])
except:
    print "Usage: roll+streak.py mesh_file roll.dat streak.dat Rey"
    sys.exit(1)

# -- Get the mesh structure information from first line of mesh_file.

mfile = open (mfilename, 'r')
nrnsnznel = mfile.readline()
if 'NR NS NZ NEL' not in nrnsnznel:
    exit('mesh_file does not have standard header')
else:
    numbers = nrnsnznel.split()
    npoints = int(numbers[0])*int(numbers[1])*int(numbers[3])
    mfile.close()

# -- Print a generic semtex field file header with this structure inserted.

print """roll+streak               Session
Tue Jan 01 00:00:00 1974  Created
""",
print '%5s%5s%5s%5s    ' % (numbers[0], numbers[1], 1, numbers[3]),
print """ Nr, Ns, Nz, Elements"""
print"""0                         Step
0                         Time
0.001                     Time step
1                         Kinvis
1                         Beta
uvwp                      Fields written
ASCII                     Format
""",

# -- Now deal with roll and streak data.

rfile = open (rfilename, 'r')
sfile = open (sfilename, 'r');

if npoints != int(rfile.readline().split()[0]):
   exit('number of points in roll file does not match declaration in mesh file')
if npoints != int(sfile.readline().split()[0]):
   exit('number of points in streak file doesnt match declaration in mesh file')

while 1:
    rline = rfile.readline();
    sline = sfile.readline()
    if not rline: break;
    rolldat = rline.split();
    stredat = sline.split();
    print str(scale*float(rolldat[2])),\
          str(scale*float(rolldat[3])),\
          stredat[2],\
          str(scale*float(rolldat[4]))

rfile.close()
sfile.close()







