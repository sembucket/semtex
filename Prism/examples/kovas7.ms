#!/bin/env mscope
#
# Kovasznay flow in the y--z plane has the exact solution:
#
#       u = 0
#       v = 1 - exp(lambda*y)*cos(2*PI*z)
#       w = lambda/(2*PI)*exp(lambda*y)*sin(2*PI*z)
#       p = (1 - exp(2*lambda*y))/2 + const.
#
# where lambda = Re/2 - sqrt(0.25*Re*Re + 4*PI*PI).
#
# We supply the exact essential BCs on (xmin,xmax), periodic on (ymin,ymax).
#
# $Id$

domain kovas7

param Re     = 40.
param DT     = 0.008
param TIME   = 1.
param LZ     = 1.
param EQTYPE = 1
param TORDER = 3
param NORDER = 8
param NZ     = 4
param lambda = Re/2-sqrt(Re^2/4+4*PI^2)
param a      = lambda

# Build the mesh 

element <-0.5, -0.5, 0.5, 0.5>
  bcond <3, X>
  bcond <0, v, 0, 1-exp(a*y)*cos(2*PI*z), a*exp(a*y)*sin(2*PI*z)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <3, X>
  bcond <2, v, 0, 1-exp(a*y)*cos(2*PI*z), a*exp(a*y)*sin(2*PI*z)/(2*PI)>

element <0, -0.5, 1, 0.5>
  bcond <1, X>
  bcond <0, v, 0, 1-exp(a*y)*cos(2*PI*z), a*exp(a*y)*sin(2*PI*z)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, X>
  bcond <2, v, 0, 1-exp(a*y)*cos(2*PI*z), a*exp(a*y)*sin(2*PI*z)/(2*PI)>

connect

# Define initial conditions

icond <Restart, kovas7.rst>

set u = 0
set v = 1-exp(a*y)*cos(2*PI*z)
set w = a*exp(a*y)*sin(2*PI*z)/(2*PI)

flowout kovas7.rst
meshout kovas7.rea

# NB: connectivity info must be generated "by hand" for the parallel code,
# because it's unsafe to rely on the system() call.  Also, the following
# assumes we're launching an MPI appliation (only system currently supported).

! mason -o kovas7.mor kovas7.rea
! mpirun -np 2 prism3d_mp -v0 kovas7

flowin kovas7.fld

enorm u 0
enorm v 1-exp(a*y)*cos(2*PI*z)
enorm w a*exp(a*y)*sin(2*PI*z)/(2*PI)

! rm kovas7.rea kovas7.fld kovas7.rst kovas7.mor


