#!/bin/env mscope
#
# Taylor flow in the x--z plane has the exact solution:
#
#       u = -cos(PI*x)*sin(PI*z)*exp(-2.0*PI*PI*KINVIS*t)
#       v =  0
#       w =  sin(PI*x)*cos(PI*z)*exp(-2.0*PI*PI*KINVIS*t)
#       p = -0.25*(cos(2.0*PI*x)+cos(2.0*PI*z))*exp(-4.0*PI*PI*KINVIS*t)
#
# Use periodic BCs.
#
# $Id$

domain taylor3

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

icond <Restart, taylor3.rst>

set u = -cos(PI*x)*sin(PI*z)
set v =  0
set w =  sin(PI*x)*cos(PI*z)

flowout taylor3.rst
meshout taylor3.rea

! prism3d -v0 taylor3 > /dev/null

flowin taylor3.fld

enorm u -cos(PI*x)*sin(PI*z)*exp(-2.0*PI*PI*TIME/Re)
enorm v  0
enorm w  sin(PI*x)*cos(PI*z)*exp(-2.0*PI*PI*TIME/Re)

! rm taylor3.rea taylor3.rst taylor3.fld
