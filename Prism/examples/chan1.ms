#
# High-resolution channel mesh
#
# $Id$

domain chan1

param Re      = 13800.
param DT      = 0.0025
param TIME    = 25.
param NORDER  = 9
param TORDER  = 2
param NZ      = 64
param IO_MEA  = 10
param IO_STAT = 10
param IO_CFL  = 50
param EQTYPE  = 1

# Mesh size in the (x,z)-plane

param XL = 2*PI
param LZ = PI

# Constant forcing in the x-direction

force <2.23e-3, 0, 0>

# Specify how to get initial conditions

icond <Restart, g3.rst>

# History points

history <uvwp, 0, 0, 0>

# Build the mesh

define yo = -1
define dy = 0.015
element <-XL/2, $yo, XL/6, $dy>
  bcond <0, W>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
     bcond <0, W>
}
bcond <1, X>

define yo = -0.985
define dy =  0.035
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = -0.95
define dy =  0.1
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = -0.85
define dy =  0.25
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = -0.6
define dy =  0.6
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = 0
define dy = 0.6
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = 0.6
define dy = 0.25
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = 0.85
define dy = 0.1
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = 0.95
define dy = 0.035
element <-XL/2, $yo, XL/6, $dy>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
}
bcond <1, X>

define yo = 0.985
define dy = 0.015
element <-XL/2, $yo, XL/6, $dy>
  bcond <2, W>
  bcond <3, X>
do j = 2, 6 {
   element <+, $yo, XL/6, $dy>
     bcond <2, W>
}
bcond <1, X>

connect

