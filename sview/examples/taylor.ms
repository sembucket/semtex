#!/bin/env mscope
#
# Visualize Taylor flow using sview
#
# Note: Run this script in "batch" mode using 'mscope -b', unless you want
#       to have mscope's display window as well as sview's.
# 
# $Id$

domain taylor

param NORDER = 11
param NZ     = 8

element <0, 0, 2, 2>
  bcond <0, Y>
  bcond <2, Y>
  bcond <1, X>
  bcond <3, X>

connect ; refine -1 ; connect

set u = -cos(PI*x)*sin(PI*z)
set v =  0
set w =  sin(PI*x)*cos(PI*z)

flowout taylor.fld
meshout taylor.rea

! meshgen taylor.rea > taylor.mesh
! ../sview -s taylor.sv taylor.mesh taylor.fld
#! rm taylor.rea taylor.mesh taylor.fld

