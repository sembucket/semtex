# -- An example SM file to produce waves

set n  = 0,32,1
set x  = 2*PI*n/32
set cx = cos(x)
set sx = sin(x)

limits x cx
box
connect x cx
connect x sx

print waves.dat { x cx sx }
