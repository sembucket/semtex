#!/bin/env mscope
#
# Channel Flow w/specified flowrate has the exact solution:
#
#     u = 1-y*y
#     v = 0
#     p = const.
#
# Flowrate is maintained at Q = 4/3.
#
# Boundary conditions are periodic in the x-direction.
#
# $Id$

domain channel2

param Re       = 100.
param DT       = 0.05
param TIME     = 50*DT
param EQTYPE   = 1
param NORDER   = 5
param NZ       = 1
param FLOWRATE = 4/3

# Build the mesh 

element <0,-1,1.5,1>
  bcond <0, W>
  bcond <3, X>
element <+,-1,1.5,1>
  bcond <0, W>
  bcond <1, X>
element <0, +,1.5,1>
  bcond <2, W>
  bcond <3, X>
element <+, 0,1.5,1>
  bcond <1, X>
  bcond <2, W>

connect

# Define initial conditions

set u = 1-y*y
set v = 0
icond <Restart, channel2.rst>
force <0,0>                    # This is required to activate FLOWRATE!

flowout channel2.rst
meshout channel2.rea

! prism2d -v0 channel2

flowin channel2

enorm u 1-y*y
enorm v 0

! rm channel2.rea channel2.rst channel2.fld
