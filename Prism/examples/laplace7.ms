#!/bin/env mscope
#
# Laplace problem on unit square, essential BC u(x,y)=sin(x)*exp(-y),
# which is also the analytical solution.  Use essential BC on upper
# edge, natural BC on other edges.
#
# $Id$

domain laplace6

param NORDER = 11

element <0, 0, 1, 1>
  bcond <0, f,  sin(x)*exp(-y)>
  bcond <1, f,  cos(x)*exp(-y)>
  bcond <2, v,  sin(x)*exp(-y)>
  bcond <3, f, -cos(x)*exp(-y)>

refine -1
connect

set f = 0
solve u f
enorm u sin(x)*exp(-y)


