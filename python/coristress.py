#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# WORKING!!!
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

if __name__ == "__main__":
    """
    """
    
    parser = argparse.ArgumentParser(description="coristress")
    parser.add_argument("alpha", type=float, help="tilt angle, deg")
    parser.add_argument("omega", type=float, help="omega")
    #parser.add_argument("-m", "--minalise", action="store_true", help="normalise fields individually using minimum value")
    parser.add_argument('infile', metavar='file', type=str, nargs='+',
                   help='the file')
    args = parser.parse_args()

    ff = fieldfile.Fieldfile(args.infile[0], "r")
    data = ff.read()
    #elmt_wise = data.reshape((ff.nr, ff.ns, ff.nel, ff.nz, ff.nflds))

    data_dict = {} 
    for i, field in enumerate(ff.fields):
        data_dict[field] = data[i]

    # -- expect GHIJKL
    #    subtract product of mean
    alpha = np.deg2rad(args.alpha)
    omega = args.omega
    Omega1 = omega
    Omega2 = (1.0-omega)/np.cos(alpha)
    mean_OMEGA_x = Omega1 + Omega2 * np.cos(alpha)
    mean_OMEGA = np.array([mean_OMEGA_x, 0., 0.])

    # -- subtract product of means, per plane

    # -- assemble time averaged Omega vector in cylindrical coordinates
    alocal_x = mean_OMEGA[0]
    alocal_y = np.zeros_like(data_dict['G']).reshape((ff.nz, ff.nel, ff.ns, ff.nr))
    alocal_z = np.zeros_like(data_dict['G']).reshape((ff.nz, ff.nel, ff.ns, ff.nhr))

    for k in range(ff.nz):
        phi = k * 2. * np.pi / ff.nz
        #sys.stderr.write('th %g\n' % phi)
        alocal_y[k,:,:,:] =  np.cos(phi) * mean_OMEGA[1] + np.sin(phi) * mean_OMEGA[2]
        alocal_z[k,:,:,:] = -np.sin(phi) * mean_OMEGA[1] + np.cos(phi) * mean_OMEGA[2]

    alocal_y = alocal_y.flatten()
    alocal_z = alocal_z.flatten()

    data_dict['G'] -= alocal_y * data_dict['w']
    data_dict['H'] -= alocal_z * data_dict['u']
    data_dict['I'] -= alocal_x * data_dict['v']

    data_dict['G'] += alocal_z * data_dict['v']
    data_dict['H'] += alocal_x * data_dict['w']
    data_dict['I'] += alocal_y * data_dict['u']

    out_u = -2. * data_dict['G']
    out_v = -2. * data_dict['H']
    out_w = -2. * data_dict['I']

    if 0:
        sys.stderr.write("u %e, %e\n" % (data_dict['u'].min(), data_dict['u'].max()))
        sys.stderr.write("v %e, %e\n" % (data_dict['v'].min(), data_dict['v'].max()))
        sys.stderr.write("w %e, %e\n" % (data_dict['w'].min(), data_dict['w'].max()))

        sys.stderr.write("G %e, %e\n" % (data_dict['G'].min(), data_dict['G'].max()))
        sys.stderr.write("H %e, %e\n" % (data_dict['H'].min(), data_dict['H'].max()))
        sys.stderr.write("I %e, %e\n" % (data_dict['I'].min(), data_dict['I'].max()))
        sys.stderr.write("ou %e, %e\n" % (out_u.min(), out_u.max()))
        sys.stderr.write("ov %e, %e\n" % (out_v.min(), out_v.max()))
        sys.stderr.write("ow %e, %e\n" % (out_w.min(), out_w.max()))

    # -- write output
    hdr = copy.copy(ff.hdr)
    hdr.fields = 'uvw'
    ff_out = fieldfile.Fieldfile('/dev/stdout', 'w', hdr)
    #ff_out.write(np.vstack(data_dict[f] for f in data_out_list))
    ff_out.write(np.vstack([out_u, out_v, out_w]))
    
    
    #    print data
    #    - of nr * ns nodes
    #    - of 11th element
    #    - of first z-plane
    #    - of second field (in this case, v)
#    print elmt_wise[:,:,10,0,1]
