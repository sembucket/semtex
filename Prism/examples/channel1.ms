#!/bin/env mscope
#
# Channel flow w/specified pressure gradient has the exact solution:
#
#         u = 1-y*y
#         v = 0
#         p = const. 
#
# Boundary conditions are periodic in the x-direction.
#
# $Id$

domain channel1

param Re     = 100.
param DT     = 0.05
param TIME   = 50*DT
param EQTYPE = 1
param NORDER = 5
param NZ     = 1

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

# Define initial conditions and forcing

set u = 1-y*y
set v = 0
icond <Restart, channel1.rst>
force <2/Re, 0>

flowout channel1.rst
meshout channel1.rea

! prism2d -v0 channel1 > /dev/null

flowin channel1

enorm u 1-y*y
enorm v 0

! rm channel1.rea channel1.rst channel1.fld

