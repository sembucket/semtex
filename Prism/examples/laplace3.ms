#!/bin/env mscope
#
# Laplace problem on unit square, with natural BC du/dn = sinPix on
# top edge.
# 
# Analytical solution is u(x,y) = sin(PIx)*sinh(PIy)/(PI*cosh(PI)).
#
# $Id$

domain laplace3

param NORDER = 12

element <0, 0, 1, 1>
  bcond <0, v, 0.>
  bcond <1, v, 0.>
  bcond <2, f, sin(PI*x)>
  bcond <3, v, 0.>

connect
  
set f = 0
solve u f
enorm u sin(PI*x)*sinh(PI*y)/(PI*cosh(PI))
