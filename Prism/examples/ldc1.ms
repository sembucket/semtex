#
# Adaptive simulation of flow in a lid-driven cavity
#
# This script computes the LDC flow using an adaptive procedure based
# on the Legendre spectrum of the (uv)-components of the velocity field. 
#
# $Id$

define nlev = 20         # number of contour levels
define tmin = 0.5        # minimum absolute value of vorticity to plot
define tmax = 10         # maximum absolute value of vorticity to plot
define lmax = 6

domain ldc

param Re     = 250
param DT     = 0.01
param NORDER = 8
param TIME   = 10
param EQTYPE = 1

icond <Restart, ldc1.rst>

# Computational domain is the unit square 

element <0,0,1,1>
  bcond <0, W, 0>
  bcond <1, W, 0>
  bcond <3, W, 0>
  bcond <2, v, 1, 0>

connect

do level = 1, 1 {
   refine -1
}

set u = y
set v = 0

meshout ldc1.rea
flowout ldc1.rst

# Generate the initial flow.  The integration here is carried out for
# a relatively long time (T=10) to get rid of the transients.  After
# this we switch to a shorter integration time for the AMR.

! echo "Computing..."
! prism2d -v0 ldc1 > /dev/null
! curl -r ldc1 ldc1.fld > ldc1.vrt

flowin ldc1.vrt

erase 
grid
bc

ltype 1 ; contour t  $tmin  $tmax $nlev
ltype 0 ; contour t -$tmax -$tmin $nlev

# Now switch to a shorter integration time and begin AMR

param TIME = 0.5

do pass = 1, 5 {
   define eps = 0.001/2^$pass

   spectrum u $eps $lmax 1
   spectrum v $eps $lmax 1
   grid

   meshout ldc1.rea
   flowout ldc1.rst

   ! echo "Computing..."
   ! prism2d -v0 ldc1 > /dev/null
   ! curl -r ldc1 ldc1.fld > ldc1.vrt

   flowin ldc1.vrt

   erase
   grid
   bc

   ltype 1 ; contour t  $tmin  $tmax $nlev
   ltype 0 ; contour t -$tmax -$tmin $nlev
}

# Cleanup

! /bin/rm ldc1.rea ldc1.fld ldc1.rst ldc1.vrt
