#!/bin/env mscope
#
# Taylor-Green vortex
#
# This script computes the evolution of the Taylor-Green vortex.
#
# $Id$

domain tg

param Re     = 3000
param NORDER = 16
param NZ     = 8
param DT     = 0.01
param TIME   = 0.1
param EQTYPE = 1
param LZ     = 2*PI
param IO_MEA = 1

# Define the mesh

element <0, 0, 2*PI, 2*PI>
  bcond <0, Y>
  bcond <2, Y>
  bcond <1, X>
  bcond <3, X>

connect
refine -1

# Initial conditions

icond <Given, sin(x)*cos(y)*cos(z), -cos(x)*sin(y)*cos(z), 0 >
force <0>

# Run the simulation

meshout tg.rea
! prism3d -v0 tg.rea
! cat tg.mea

# Clean up

! rm tg.rea tg.fld tg.mea




