#!/bin/env mscope
#
# Laplace problem on unit square, essential BC u(x,y)=sin(x)*exp(-y),
# which is also the analytical solution.
#
# Note that the mesh has one curved side. 
#
# $Id$

domain laplace5

param NORDER = 11

element <0, 0, 1, 1>
  bcond <0, v, sin(x)*exp(-y)>
  bcond <1, v, sin(x)*exp(-y)>
  bcond <2, v, sin(x)*exp(-y)>
  bcond <3, v, sin(x)*exp(-y)>

refine -1

set f = 0
solve u f
enorm u sin(x)*exp(-y)


