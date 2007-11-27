lweight 3
ltype 0
set alpha = 0,2,0.01
set Amax  = 1 - 1.177*alpha + 0.351*alpha*alpha

limits 0 2 0 1
location 5000 31000 5000 21000
expand 1.5
ticksize 0.2 1 0.1 0.5
box
expand 1.8 
xlabel \alpha
ylabel {\i A}^{\ast}
#connect alpha Amax

data scanned.txt
read {a 1 b 2 }
ptype 20 3
expand 1.5
points a b

set x = { 0  0.1    0.25   0.5    0.75   1      1.25    1.5     2 }
set y = { 1.0044 0.9004 0.7542 0.4927 0.2718 0.1490 0.1032 0.0821 0.0618 }
set sx = 0,2,0.01
ptype 4 1
expand 4
spline x y sx sy
connect sx sy
#points x y 

lweight 2
relocate 0.9514 0 draw 0.9514 0.25
expand 1.2
relocate 0.85 0.30 label 1.0514^{-1}
ltype 1
relocate 0 1 draw 0.9514 0




