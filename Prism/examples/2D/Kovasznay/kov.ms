#!/bin/env mscope
#
# Kovasznay Flow
#

domain kov

param Re     = 40.
param DT     = 0.001
param TIME   = 50*DT
param EQTYPE = 1
param TORDER = 3
param NORDER = 8
param lambda = Re/2-sqrt(Re^2/4+4*PI^2)  # flow parameter
param a      = lambda                    # shorthand

# Build the mesh 

element <-0.5, -0.5, 0.5, 0.5>
  bcond <0, Y>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <-0.5, +, 0.5, 0.5>
  bcond <2, Y>
  bcond <3, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>

element <0, -0.5, 1, 0.5>
  bcond <0, Y>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
element <0, +, 1, 0.5>
  bcond <1, v, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>
  bcond <2, Y>

connect

# Define initial conditions

icond <Given, 1-exp(a*x)*cos(2*PI*y), a*exp(a*x)*sin(2*PI*y)/(2*PI)>

# Save to the file kov.rea

meshout kov.rea

# Build program and run convergence study

! gmake > /dev/null

define pmin = 4
define pmax = 12

! echo "Convergence study [NORDER $pmin to $pmax]:"
do n = $pmin, $pmax {
   ! ./kov2d -n $n kov | tail -1
}
! gmake empty > /dev/null
! /bin/rm kov.rea
