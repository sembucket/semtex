#!/bin/env mscope
#
# Poisson problem on the unit square, with homogeneous essential BCs
# and forcing f = -sin(PI*x)*sin(PI*y).
#
# Analytical solution is u(x,y) = -1.0 / (2*PI^2) * sin(PI*x)*sin(PI*y)
#
# $Id$

domain poisson1

param NORDER = 11

element <0, 0, 1, 1>
  bcond <0, v, 0.>
  bcond <1, v, 0.>
  bcond <2, v, 0.>
  bcond <3, v, 0.>
  
refine -1

set f = -sin(PI*x)*sin(PI*y)
solve u f
enorm u -1.0/(2*PI*PI)*sin(PI*x)*sin(PI*y)


