#!/bin/env mscope
#
# Couette Flow
# 
# This flow depends on three parameters: the velocity of the upper wall
# (Uo), the  pressure gradient (dP), and the Reynolds number (Re).
# The reference scales for the problem are Uo and the channel height
# (h=1).  The Reynolds number is defined in terms of these param-
# meters as Re=Uo*h/\nu.  The leaves dP and Re as the only real
# "control parameters" in the problem, with the non-dimensional upper-
# wall velocity equal to 1.
#
# When you run this script, it creates a plot of the steady-state
# velocity profile.  The exact and computed profiles are overlayed,
# but hopefully you will only see one curve!
#
# $Id$

domain couette

# The two adjustable parameters:

param dP = 0.05  # Pressure gradient, dP > 0 (adverse), dP < 0 (favorable) 
param Re = 100.  # Reynolds number

# Other parameters...

param Uo     = 1.
param DT     = 0.05
param TIME   = 10*DT
param EQTYPE = 1
param NORDER = 5

force <2, -dP, 0>   # Note: body force = -(dp/dx)

# Build the mesh 

element <0, 0,1.5,1>
  bcond <0, W>
  bcond <2, V, Uo, 0.>
  bcond <3, X>
element <+, 0,1.5,1>
  bcond <0, W>
  bcond <1, X>
  bcond <2, V, Uo, 0.>

connect

# Compute the exact solution and use it as the initial condition

set u = y*Uo-0.5*dP*Re*y*(1-y)
set v = 0
flowout couette.rst
icond <Restart, couette.rst>

# Run the simulation and plot a profile

option autoProfile on
ltype 0 ; profile u 0 0 0 1

meshout couette.rea
! prism2d -v0 couette 
flowin  couette.fld

option autoProfile off
ltype 1 ; profile u 0 0 0 1

enorm u y*Uo-0.5*dP*Re*y*(1-y)

# Clean up

! rm couette.rea couette.fld couette.rst

