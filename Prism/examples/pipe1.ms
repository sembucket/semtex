#!/bin/env mscope
#
# Pipe flow w/specified pressure gradient has the exact solution:
#
#    u = 0
#    v = 0
#    w = 1-(x^2+y^2)
#
# Constant forcing of f=(0,0,4/Re) is applied to maintain steady flow.
# Boundary conditions are solid walls around the perimeter of the pipe 
# and periodic in the streamwise z-direction.
#
# $Id$

domain pipe1

param Re     = 100.
param DT     = 0.05
param TIME   = 50*DT
param EQTYPE = 1
param NORDER = 8
param NZ     = 4

# Build the mesh 

define r = 0.5
node <n0, $r*cos(  PI/4), $r*sin(  PI/4)>
node <n1, $r*cos(3*PI/4), $r*sin(3*PI/4)>
node <n2, $r*cos(5*PI/4), $r*sin(5*PI/4)>
node <n3, $r*cos(7*PI/4), $r*sin(7*PI/4)>

define r = 1.0
node <n4, $r*cos(  PI/4), $r*sin(  PI/4)>
node <n5, $r*cos(3*PI/4), $r*sin(3*PI/4)>
node <n6, $r*cos(5*PI/4), $r*sin(5*PI/4)>
node <n7, $r*cos(7*PI/4), $r*sin(7*PI/4)>

# center element
element <n2, n3, n0, n1>

# arc segments
element <n1, n0, n4, n5>
  curve <2, arc, 1.0>
  bcond <2, W>
element <n6, n2, n1, n5>
  curve <3, arc, 1.0>
  bcond <3, W>
element <n6, n7, n3, n2>
  curve <0, arc, 1.0>
  bcond <0, W>
element <n3, n7, n4, n0>
  curve <1, arc, 1.0>
  bcond <1, W>

connect

# Initial conditions

force <0, 0, 4/Re>
icond <Restart, pipe1.rst>

set u = 0
set v = 0
set w = 1-(x*x+y*y)

flowout pipe1.rst
meshout pipe1.rea

! prism3d pipe1 > /dev/null

flowin pipe1.fld

enorm u 0
enorm v 0
enorm w 1-(x*x+y*y)

! rm pipe1.rea pipe1.rst pipe1.fld
