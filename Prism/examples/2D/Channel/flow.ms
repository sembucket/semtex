#!/bin/env mscope
#
# Channel Flow
#
# This is similar to the test case "chan.ms" except that the flow is driven by
# specifying a flowrate rather than a body force (pressure gradient).  This is
# controlled by setting the parameter FLOWRATE to the desired value, as 
# defined by:
#
#               FLOWRATE = -\int_B u*n ds,
# 
# where B is the periodic boundary segment of the domain, u[x,y] is the velo-
# city field, and n is the unit outward normal.
#
# This script sets exact initial conditions so the velocity field should not
# change in time.
#
# $Id$

domain chan

param Re       = 100.
param DT       = 0.05
param TIME     = 10.
param EQTYPE   = 1
param NORDER   = 5
param FLOWRATE = 4./3.
param IO_MEA   = 1

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

icond <Given, 1-y*y, 0>
force <2, 0, 0>

# Define probes

probe <uvp, 0., 0., 0>
probe <uvp, 1., 0., 0>
probe <uvp, 2., 0., 0>

# Run the simulation and plot a profile across the channel

meshout chan.rea
! prism2d -v0 chan 
flowin chan.fld
erase
profile u 0 -1 0 2

! rm chan.rea chan.his chan.mea chan.fld
