#!/bin/env mscope
#
# Kovasznay flow has the exact solution:
#
#       u = 1 - exp(lambda*x)*cos(2*PI*y)
#       v = lambda/(2*PI)*exp(lambda*x)*sin(2*PI*y)
#       p = (1 - exp(2*lambda*x))/2 + const.
#
# where lambda = Re/2 - sqrt(0.25*Re*Re + 4*PI*PI).
#
# We supply the exact essential BCs on velocity.
#
# $Id$

domain kovas1

param Re     = 40.
param DT     = 0.008
param TIME   = 1.
param EQTYPE = 1
param TORDER = 3
param NORDER = 8
param NZ     = 1
param lambda = Re/2-sqrt(Re^2/4+4*PI^2)
param a      = lambda

# Build the mesh 

element <-0.5, -0.5, 0.5, 0.5>
  bcond <0, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <2, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>

element <0, -0.5, 1, 0.5>
  bcond <0, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
  bcond <2, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>

connect

# Define initial conditions

set u = 1-exp(a*x)*cos(2*PI*y)
set v = a*exp(a*x)*sin(2*PI*y)/(2*PI)

icond <Restart, kovas1.rst>
flowout kovas1.rst
meshout kovas1.rea

! prism2d -v0 kovas1

flowin kovas1.fld

enorm u 1-exp(a*x)*cos(2*PI*y)
enorm v a*exp(a*x)*sin(2*PI*y)/(2*PI)

! rm kovas1.rea kovas1.fld kovas1.rst
