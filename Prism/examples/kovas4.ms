#!/bin/env mscope
#
# Kovasznay flow in the x--z plane has the exact solution:
#
#       u = 1 - exp(lambda*x)*cos(2*PI*z)
#       v = 0
#       w = lambda/(2*PI)*exp(lambda*x)*sin(2*PI*z)
#       p = (1 - exp(2*lambda*x))/2
#
# where lambda = Re/2 - sqrt(0.25*Re*Re + 4*PI*PI).
#
# We supply the exact essential BCs on (xmin,xmax), periodic on (ymin,ymax)
#
# $Id$

domain kovas4

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
  bcond <0, Y>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*z), 0, a*exp(a*x)*sin(2*PI*z)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <2, Y>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*z), 0, a*exp(a*x)*sin(2*PI*z)/(2*PI)>

element <0, -0.5, 1, 0.5>
  bcond <0, Y>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*z), 0, a*exp(a*x)*sin(2*PI*z)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*z), 0, a*exp(a*x)*sin(2*PI*z)/(2*PI)>
  bcond <2, Y>

connect

# Define initial conditions

icond <Restart, kovas4.rst>

set u = 1-exp(a*x)*cos(2*PI*z)
set v = 0
set w = a*exp(a*x)*sin(2*PI*z)/(2*PI)

flowout kovas4.rst
meshout kovas4.rea

! prism3d -v0 kovas4

flowin kovas4.fld

enorm u 1-exp(a*x)*cos(2*PI*z)
enorm v 0
enorm w a*exp(a*x)*sin(2*PI*z)/(2*PI)

! rm kovas4.rea kovas4.fld kovas4.rst
