#!/bin/env mscope
#
# Channel Flow
#
# This script runs a channel flow simulation at Re=100, starting from default
# initial conditions.  The flow is driven by prescribing a body force (pres-
# sure gradient) in the x-direction.
#
# In the plotting window you should see a sequence of velocity pro-
# files showing the developing flow, eventually reaching the steady-
# state solution u(y)=y*(1-y).
#
# $Id$

# Number of sequences to run.  Feel free to change this, but if you
# lower it the flow will probably not reach the steady-state
# solution. 

define nmax  = 50

domain chan

param Re     = 100.
param DT     = 0.05
param TIME   = 5.
param EQTYPE = 1
param NORDER = 5

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

set u = 0
set v = 0
icond <Restart, chan.rst>
force <2, 12/Re, 0>

# Run the simulation and plot profiles

meshout chan.rea

option autoProfile off
erase
limits 0 2 -.1 1.1
box
ltype 2
relocate 0 1 ; draw 2 1
ltype 0

do n = 1, $nmax {
   flowout chan.rst
   ! prism2d -v0 chan
   flowin  chan.fld
   profile u 0 -1 0 2
}

# Clean up

! rm chan.rea chan.fld chan.rst

