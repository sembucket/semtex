#!/bin/env mscope
#
# Plane Poiseulle flow on a nonconforming mesh
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

domain channel3

param Re       = 100.
param DT       = 0.01
param TIME     = 50*DT
param EQTYPE   = 0
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

refine -1
refine  4  5  8  9
refine 14 15 18 19

# Define initial conditions

set u = 1-y*y
set v = 0
icond <Restart, channel3.rst>
force <0,0>                    # This is required to activate FLOWRATE!

flowout channel3.rst
meshout channel3.rea

! prism2d -v0 channel3

flowin channel3

enorm u 1-y*y
enorm v 0

! rm channel3.rea channel3.rst channel3.fld
