#!/usr/bin/env python
#
# Usage: packwave.py mesh_file wave.dat
# 
# Take a wave.dat file in Spen's format and rearrange it so it has the
# same format as a dog complex eigenmode (ASCII).  The input data is
#
# u.Re I.Im v.Re v.Im w.Re w.Im p.Re p.Im
#
# at each mesh point, whereas we want
#
# u.Re v.Re w.Re p.Re followed by
# u.Im v.Im w.Im p.Im
#
# The (meshpr-type) mesh file is used to generate structure
# information that gets used to make a semtex-style header.
# 
# $Id$
##############################################################################

import sys, os

try:
    mfilename = sys.argv[1]
    wfilename = sys.argv[2]
except:
    print "Usage: packwave.py mesh_file wave.dat"
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
print '%5s%5s%5s%5s    ' % (numbers[0], numbers[1], 2, numbers[3]),
print """ Nr, Ns, Nz, Elements"""
print"""0                         Step
0                         Time
0.001                     Time step
1                         Kinvis
1                         Beta
uvwp                      Fields written
ASCII                     Format
""",

# -- Now deal with wave data.

wfile = open (wfilename, 'r')

if npoints != int(wfile.readline().split()[0]):
   exit('number of points in wave file does not match declaration in mesh file')

# -- First all the real data.

while 1:
    wline = wfile.readline();
    if not wline: break;
    wavedat = wline.split();
    print wavedat[2], wavedat[4], wavedat[6], wavedat[8]

# -- Rewind and do all the imaginary data.

wfile.seek(0);
wfile.readline();

while 1:
    wline = wfile.readline();
    if not wline: break;
    wavedat = wline.split();
    print wavedat[3], wavedat[5], wavedat[7], wavedat[9]

wfile.close()








