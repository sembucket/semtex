#!/bin/env mscope
#
# Laplace problem on unit square, with homogeneous essential BCs
# on lower, left- & right-hand edges and u(x) = sin(2PIx) on upper edge.
#
# Analytical solution is c(x,y) = sin(2PIx)*sinh(2PIy)/sinh(2PI).

domain laplace1

param NORDER = 11

element <0, 0, 1, 1>
  bcond <0, v, 0.>
  bcond <1, v, 0.>
  bcond <2, v, sin(2*PI*x)>
  bcond <3, v, 0.>
  
refine -1

set f = 0
solve u f
enorm u sin(2*PI*x)*sinh(2*PI*y)/sinh(2*PI)


