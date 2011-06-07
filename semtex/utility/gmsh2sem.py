#!/usr/bin/env python
#
# Read a gmsh (http://www.geuz.org/gmsh/) mesh file (.msh extension
# typical but we do not check), assumed to supply an all-quad mesh,
# and from this, produce a semtex session file skeleton, containing
# NODES and ELEMENTS sections. Wrap this with other standard semtex
# feml file sections in order to generate a valid semtex session file.
#
# Print all this to standard output.
#
# Usage: gmsh2sem.py gmshfile.msh
# Usage: python gmsh2sem.py gmshfile.msh
#
# NB: use Mesh.SubdivisionAlgorithm=1; in gmsh in order to generate quads.
# Standard semtex utility meshpr must be found in your shell's PATH variable.
# 
# $Id$
##############################################################################

import sys, os

try:
    infilename = sys.argv[1]
except:
    print "Usage: gmesh2sem infile.msh"; sys.exit(1)

ifile = open (infilename, 'r')
t1    = open (".t1", 'w')


# -- We will assume that the file is a valid gmsh mesh file with quad meshes.
#    First print a generic semtex header.

s = """# Skeleton semtex session file generated from gmsh .msh file

<TOKENS>
        N_P    = 5
        N_STEP = 20
        D_T    = 0.01
        KINVIS = 1.
</TOKENS>

<FIELDS NUMBER=3>
        u v p
</FIELDS>

<GROUPS NUMBER=1>
        1       w       wall
</GROUPS>

<BCS NUMBER=1>
        1       w       3
                <D> u = 0. </D>
                <D> v = 0. </D>
                <H> p      </H>
</BCS>\n
"""

t1.write(s)

# -- Extract NODES and ELEMENTS from given gmsh file.

while 1:
    line = ifile.readline()
    if '$Nodes' in line:
        line = ifile.readline()
        Nnode = int (line)
        # -- Nodes are written to t1 on-the-fly
        s = '<NODES NUMBER=' + repr(Nnode) + '>\n'
        t1.write(s)
        for i in range (Nnode):
            line = ifile.readline()
            t1.write(line),
        t1.write('</NODES>\n\n')
    elif '$Elements' in line:
        # -- Since we don't know the number of valid (i.e. quad) elements yet,
        #    we store them in a list first and write to t1 later on.
        Elist = []
        line = ifile.readline()
        Nel = int (line)
        Nquads = 0
        for i in range (Nel):
            line = ifile.readline()
            words = line.split()
            # -- A gmsh file may contain non-quad elements. Skip those.
            if int(words[1]) != 3: continue
            Nquads += 1
            # -- Quad element lines have 10 numbers.
            #    The first one is a tag, and the last 4 are node numbers.
            s = '%i <Q> %s %s %s %s </Q>\n' % \
                (Nquads, words[6], words[7], words[8], words[9])
            Elist.append(s)
        Elist.append('</ELEMENTS>\n\n')

        break

ifile.close()

# -- Add number of (quad) elements and write element list to file
Elist.insert(0, '<ELEMENTS NUMBER=%i>\n' % Nquads)
for s in Elist:
    t1.write(s)
t1.close()

# -- Generate a SURFACES section that contains all the unfilled
#    surface valencies (those element edges that do not mate to
#    another element). You can then edit this (by hand, or other
#    means).

os.system("meshpr -s .t1 > .t2")

t1 = open (".t1", 'r')
t2 = open (".t2", 'r')

# -- Print up all of t1 and t2

for line in t1.readlines(): print line,
for line in t2.readlines(): print line,

t1.close();       t2.close()
os.remove('.t1'); os.remove('.t2')
