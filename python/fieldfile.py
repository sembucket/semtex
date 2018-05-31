#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# python classes for accessing semtex fieldfile
# by Thomas Albrecht
#
# develop at /home/albrecht/tud/Semtex/semtex-work/utils

import string
import numpy as np
import sys

# ----------------------------------------------------------
class Geometry:
    def __init__(self, nr=3, ns=3, nz=1, nel=1):
        self.nr  = nr
        self.ns  = ns
        self.nz  = nz
        self.nel = nel
    def __str__(self):
        return "%i %i %i %i" % (self.nr, self.ns, self.nz, self.nel)

    def __eq__(self, other):
        is_equal = True
        msg = ""
        if self.nr != other.nr:
            msg += " nr differs"
            is_equal = False
        if self.ns != other.ns:
            msg += " ns differs"
            is_equal = False
        if self.nz != other.nz:
            msg += " nz differs"
            is_equal = False
        if self.nel != other.nel:
            msg += " nel differs"
            is_equal = False

        return is_equal, msg

class Header:
#    def __init__(self, nr=3, ns=3, nz=1, nel=1, geometry=None, fields="", format="binary"):
    def __init__(self, geometry, fields="", format="binary"):
        self.session = ""
        self.created = ""
        #if not geometry: geometry = Geometry(nr, ns, nz, nel)
        self.geometry = geometry
        self.step    = 0
        self.time    = 0.
        self.dt      = 0.
        self.kinvis  = 1.
        self.beta    = 1.
        self.fields  = fields
        self.format  = format


    def read(self, f):
        hdr = []
        for i in range(0,9):
            hdr.append(string.split(f.readline()))
        hdr.append(f.readline())

        self.session = hdr[0][0]
        self.created = hdr[1][0]
        self.geometry = Geometry(int(hdr[2][0]), int(hdr[2][1]), \
                                 int(hdr[2][2]), int(hdr[2][3]))
        self.step    =   int(hdr[3][0])
        self.time    = float(hdr[4][0])
        self.dt      = float(hdr[5][0])
        self.kinvis  = float(hdr[6][0])
        self.beta    = float(hdr[7][0])
        self.fields  =       hdr[8][0]
        self.format  =       hdr[9][:-1]

    def write(self, f):
        f.write(str(self))

    def __str__(self):
        out = []
        out.append("%-25s Session"        % self.session)
        out.append("%-25s Created"        % self.created)
        out.append("%-4i %-4i %-4i %-6i     Nr, Ns, Nz, Elements" % \
            (self.geometry.nr, self.geometry.ns, self.geometry.nz, self.geometry.nel))
        out.append("%-25i Step"           % self.step)
        out.append("%-25g Time"           % self.time)
        out.append("%-25g Time step"      % self.dt)
        out.append("%-25g Kinvis"         % self.kinvis)
        out.append("%-25g Beta"           % self.beta)
        out.append("%-25s Fields written" % self.fields)
        out.append("%-25s Format"         % self.format)
        return string.join(out, '\n') + '\n'


# ------------------------------------------------------------------------------
class Fieldfile:
    def __init__(self, fname, state, header=None):

        if state == "r":
            self.f = open(fname,  "r")
            self.hdr = Header(Geometry())
            self.hdr.read(self.f)

        elif state == "w":
            if not header: 
                raise ValueError('Need header when writing file')
            self.hdr = header
            self.f = open(fname,  "w")
            self.hdr.write(self.f)

        # FIXME: this design sucks. Shouldnt replicate data.
        #        Get rid of hdr?
        self.nr      = self.hdr.geometry.nr
        self.ns      = self.hdr.geometry.ns
        self.nz      = self.hdr.geometry.nz
        self.nel     = self.hdr.geometry.nel

        self.nrns    = self.hdr.geometry.nr * self.hdr.geometry.ns
        self.nxy     = self.nrns * self.hdr.geometry.nel
        self.ntot    = self.nxy  * self.hdr.geometry.nz
        self.nflds   = len(self.hdr.fields)
        self.ntotf   = self.ntot * self.nflds
        
        self.data = None

        # -- create list of field variables
        self.fields = [f for f in self.hdr.fields]

    # --------------------------------------------------------------------------
#    def read(self):
#        bin = array.array('d')
#        bin.read(self.f, 1)
#        return bin[0]

    # --------------------------------------------------------------------------
    def write(self, data, keep_open=False):
        """write field data to file"""
        if data.dtype != np.dtype('float64'):
            raise TypeError('need float64 data')
        data.tofile(self.f)
        if not keep_open:
            self.f.close()
        

    def write_fields(self, fields):
        self.write(np.vstack(the_field.flatten() for the_field in fields))

    def __getitem__(self, fieldname):
        return self.data[self.field_index(fieldname)]

    # --------------------------------------------------------------------------
    def read(self):
        "read field data from file. Float64 data expected."
#        return np.fromfile(self.f, 'd').reshape(self.nflds, self.ntot) # fails with stdin
        #return np.frombuffer(self.f, 'd', count=-1).reshape(self.nflds, self.ntot)
        buf = self.f.read()
        self.data = np.fromstring(buf, np.float64, count=-1).reshape(self.nflds, self.ntot)
        return self.data
        #return np.fromfile(self.f, 'd')
        self.close()

    def alloc(self):
        """allocate and return data storage"""
#        if not fields: fields = self.hdr.fields
        return np.zeros((self.ntotf), dtype = 'float64')

    # --------------------------------------------------------------------------
    def close(self):
        self.f.close();

    def field_index(self, needle):
        return self.fields.index(needle)

    def is_mesh_compatible(self, mesh):
        return mesh.geometry == self.hdr.geometry

# ----------------------------------------------------------
def convert():
    """convert fieldfile to ASCII, like utility/convert.C"""
    ff = Fieldfile("example.fld", "r")
    ff.hdr.write(sys.stdout)
    data = ff.read()

    for i in range(ff.ntot):
        for field in range(ff.nflds):
            print "%g" % data[i, field],
        print


def element_wise():
    """demonstrates element-wise access
       NB: Fieldfile.read() expects double precision (Float64) data, as is standard for
           semtex field files. No checks done for single precision or funny byte order.
    """
    ff = Fieldfile("example.fld", "r")
    data = ff.read()
    #elmt_wise = data.reshape((ff.nr, ff.ns, ff.nel, ff.nz, ff.nflds))
    elmt_wise = data.reshape((ff.nflds, ff.nz, ff.nel, ff.ns, ff.nr)) # works!
    #    print data
    #    - of nr * ns nodes
    #    - of 11th element
    #    - of first z-plane
    #    - of second field (in this case, v)
    print elmt_wise[1,0,10,:,:]

if __name__ == "__main__":
    element_wise()
    #convert()
    #main()
