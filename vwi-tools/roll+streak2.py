#!/usr/bin/env python
#
# Usage: roll+streak.py roll.asc streak.asc Rey
# 
# Combine u, v, p in roll file with w in streak file in order to
# generate a 2D/3C roll + streak flow.  The final argument, Rey, is
# used to scale the roll variables relative to the steak variables;
# typically this scaling is 1/Rey.
#
# In this version the two .fld files are standard semtex ASCII-format field
# files.  We check that the headers say there are the same number of data in
# each.
# 
# $Id$
##############################################################################

import sys, os

try:
    rfilename = sys.argv[1]
    sfilename = sys.argv[2]
    Rey = sys.argv[3]
    scale = 1.0/float(Rey)
except:
    print "Usage: roll+streak2.py roll.asc streak.asc Rey"
    sys.exit(1)

# -- Open the streak file, repeat header, check it's ASCII, etc.

sfile = open (sfilename, 'r')
line = sfile.readline()

print "roll+streak               Session"

line = sfile.readline()
print line,

line = sfile.readline()
print line,

numbers = line.split()
ntot = int(numbers[0])*int(numbers[1])*int(numbers[3])

line = sfile.readline()
print line,
line = sfile.readline()
print line,
line = sfile.readline()
print line,

line = sfile.readline()
print '%-25g Kinvis' % scale

line = sfile.readline()
print line,

line = sfile.readline()
format = line.split()[0]
if format != 'uvwp':
    print 'ERROR: variables in streak file must be uvwp:', format
    sys.exit(1)
print line,

line = sfile.readline()
format = line.split()[0]
if format != 'ASCII':
    print 'ERROR: Invalid format in streak file: ', format
    sys.exit(1)
print line,

# -- That's done with the streak file header, now the roll file.

rfile = open (rfilename, 'r')

line = rfile.readline()
line = rfile.readline()
line = rfile.readline()

numbers = line.split()
ntot2 = int(numbers[0])*int(numbers[1])*int(numbers[3])

if ntot != ntot2:
    print 'ERROR: number of roll and streak points differ:', ntot, ntot2
    sys.exit(1)

line = rfile.readline()
line = rfile.readline()
line = rfile.readline()
line = rfile.readline()
line = rfile.readline()
line = rfile.readline()

format = line.split()[0]
if format != 'uvp':
    print 'ERROR: variables in roll file must be uvp:', format
    sys.exit(1)

line = rfile.readline()
format = line.split()[0]
if format != 'ASCII':
    print 'ERROR: Invalid format in roll file: ', format
    sys.exit(1)
    
# -- Now deal with roll and streak data.

while 1:
    rline = rfile.readline();
    sline = sfile.readline()
    if not rline: break;
    rolldat = rline.split();
    stredat = sline.split();
    print str(scale*float(rolldat[0])),\
          str(scale*float(rolldat[1])),\
          stredat[2],\
          str(scale*float(rolldat[2]))

rfile.close()
sfile.close()
