define vT 1.92963029109e-2
define shrink 4.24990562851e-2

set dt = { \
0.02 \
0.01 \
0.005 \
0.002 \
0.001 \
0.0005 \
0.0002 \
}

set v103 = { \
2.081864e-02 \
2.005286e-02 \
1.967341e-02 \
1.944686e-02 \
1.937154e-02 \
1.933391e-02 \
1.931134e-02 \

}

set v203 = { \
1.927869e-02 \
1.929197e-02 \
1.929523e-02 \
1.929613e-02 \
1.929626e-02 \
1.929629e-02 \
1.929630e-02 \
}

set v301 = { \
4.544766e-01 \
4.541518e-01 \
4.540687e-01 \
4.540452e-01 \
4.540419e-01 \
4.540410e-01 \
4.540408e-01 \
}

set v303 = { \
1.931581e-02 \
1.930114e-02 \
1.929751e-02 \
1.929650e-02 \
1.929635e-02 \
1.929631e-02 \
1.929630e-02 \
}

set ldt = lg(dt)
set le1 = lg(abs(v103 - $vT))
set le2 = lg(abs(v203 - $vT))
#set le3 = lg(abs(v303 - $vT))
set le3 = lg(abs(v303/v301 - $shrink))

lweight 2
limits -7 -1 -9  -3
ticksize -1 0 -1 0
notation 1 1 1 1
expand 1.01
location 15000 30000 15000 30000
box
expand 1.5
xlabel \Delta \i \,t
ylabel \parallel\,decay rate\,\parallel_\infty 
ptype 4 0
points ldt le1
ptype 20 0
points ldt le2
ptype 3 0
points ldt le3
relocate $(lg(4e-3)) $(lg(abs(1.929707e-02/4.540587e-01 - $shrink))) dot
relocate $(lg(8e-3)) $(lg(abs(1.929940e-02/4.541121e-01 - $shrink))) dot
relocate $(lg(4e-2)) $(lg(abs(1.937570e-02/4.557142e-01 - $shrink))) dot


