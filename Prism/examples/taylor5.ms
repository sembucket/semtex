#!/bin/env mscope
#
# Taylor flow in the y--z plane has the exact solution
#
#       u =  0
#       v = -cos(PI*y)*sin(PI*z)*exp(-2.0*PI*PI*KINVIS*t)
#       w =  sin(PI*y)*cos(PI*z)*exp(-2.0*PI*PI*KINVIS*t)
#       p = -0.25*(cos(2.0*PI*y)+cos(2.0*PI*z))*exp(-4.0*PI*PI*KINVIS*t)
#
# Use periodic BCs.
#
# $Id$

domain taylor5

param Re     = 100.
param DT     = 0.02
param TIME   = 20*DT
param EQTYPE = 1
param TORDER = 2
param NORDER = 11
param NZ     = 8
param LZ     = 2

element <0, 0, 2, 2>
  bcond <0, Y>
  bcond <2, Y>
  bcond <1, X>
  bcond <3, X>

connect ; refine -1 ; connect

icond <Restart, taylor5.rst>

set u =  0
set v = -cos(PI*y)*sin(PI*z)
set w =  sin(PI*y)*cos(PI*z)

flowout taylor5.rst
meshout taylor5.rea

# NB: connectivity info must be generated "by hand" for the parallel code,
# because it's unsafe to rely on the system() call.  Also, the following
# assumes we're launching an MPI appliation (only system currently supported).

! mason -o taylor5.mor taylor5
! mpirun -np 4 prism3d_mp taylor5

flowin taylor5.fld

enorm u  0
enorm v -cos(PI*y)*sin(PI*z)*exp(-2.0*PI*PI*TIME/Re)
enorm w  sin(PI*y)*cos(PI*z)*exp(-2.0*PI*PI*TIME/Re)

! rm taylor5.rea taylor5.rst taylor5.fld taylor5.mor
