#!/bin/env mscope
#
# Kovasznay flow in the x--y plane has the exact solution:
#
#       u = 1 - exp(lambda*x)*cos(2*PI*y)
#       v = lambda/(2*PI)*exp(lambda*x)*sin(2*PI*y)
#       w = 0
#       p = (1 - exp(2*lambda*x))/2 + const.
#
# where lambda = Re/2 - sqrt(0.25*Re*Re + 4*PI*PI).
#
# We supply the exact essential BCs on (xmin,xmax), periodic on (ymin,ymax)
#
# $Id$

domain kovas3

param Re     = 40.
param DT     = 0.008
param TIME   = 1.
param EQTYPE = 1
param TORDER = 3
param NORDER = 8
param NZ     = 4
param lambda = Re/2-sqrt(Re^2/4+4*PI^2)
param a      = lambda

# Build the mesh 

element <-0.5, -0.5, 0.5, 0.5>
  bcond <0, Y>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI), 0>
element <-0.5, +, 0.5, 0.5>
  bcond <2, Y>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI), 0>

element <0, -0.5, 1, 0.5>
  bcond <0, Y>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI), 0>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI), 0>
  bcond <2, Y>

connect

# Define initial conditions

icond <Restart, kovas3.rst>

set u = 1-exp(a*x)*cos(2*PI*y)
set v = a*exp(a*x)*sin(2*PI*y)/(2*PI)
set w = 0

flowout kovas3.rst
meshout kovas3.rea

! prism3d -v0 kovas3

flowin kovas3.fld

enorm u 1-exp(a*x)*cos(2*PI*y)
enorm v a*exp(a*x)*sin(2*PI*y)/(2*PI)
enorm w 0

! rm kovas3.rea kovas3.fld kovas3.rst
