ticksize 0 0 -1 0
limits 1 14 -11 0
expand 1.4
notation -4 4 1 1
location 5000 31000 5000 31000
box
xlabel {\i N_p}
ylabel ||{\i u-u_e}||_\infty
relocate 5 -1
label "Offset-axis Kovasznay flow"
relocate 5 -1.5
label "Cylindrical-coordinate code"

set np = { 2 3 4 5 6 7 8 9 10 11 12 13 }
set ue = { 0.438416 0.0632847 0.0175626 0.00151925 0.000395325 1.51193e-05 2.95626e-05 2.31721e-05 1.75045e-05 1.34045e-05 1.06364e-05 8.48367e-06}

expand 2
set le = lg(ue)
ptype 4 0
points np le
relocate 3 -8
dot
relocate 3.2 -8.1
expand 1.4
label "Axis included"

set np = { 2 3 4 5 6 7 8 9 10 11 12 13 }
set ue = { 0.435102 0.0625267 0.0139614 0.00101457 0.000288583 1.37665e-05 3.39491e-06 1.15404e-07 2.25294e-08 6.12904e-10 1.48686e-10 1.48686e-10 }

expand 2
set le = lg(ue)
ptype 4 3
points np le
relocate 3 -9
dot
expand 1.4
relocate 3.2 -9.1
label "Axis excised"

