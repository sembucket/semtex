#!/bin/env mscope
#
# Taylor flow in the x--y plane has the exact solution:
#
#       u = -cos(PI*x)*sin(PI*y)*exp(-2.0*PI*PI*KINVIS*t)
#       v =  sin(PI*x)*cos(PI*y)*exp(-2.0*PI*PI*KINVIS*t)
#       w =  0
#       p = -0.25*(cos(2.0*PI*x)+cos(2.0*PI*y))*exp(-4.0*PI*PI*KINVIS*t)
#
# Use periodic BCs.
#
# $Id$

domain taylor2

param Re     = 100.
param DT     = 0.02
param TIME   = 20*DT
param EQTYPE = 1
param TORDER = 2
param NORDER = 11
param NZ     = 8
param LZ     = 2

element <0, 0, 2, 2>
  bcond <0, Y>
  bcond <2, Y>
  bcond <1, X>
  bcond <3, X>

connect ; refine -1 ; connect

icond <Restart, taylor2.rst>

set u = -cos(PI*x)*sin(PI*y)
set v =  sin(PI*x)*cos(PI*y)
set w =  0

flowout taylor2.rst
meshout taylor2.rea

! prism3d -v0 taylor2

flowin taylor2.fld

enorm u -cos(PI*x)*sin(PI*y)*exp(-2.0*PI*PI*TIME/Re)
enorm v  sin(PI*x)*cos(PI*y)*exp(-2.0*PI*PI*TIME/Re)
enorm w  0

! rm taylor2.rea taylor2.rst taylor2.fld
