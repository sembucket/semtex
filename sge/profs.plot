erase
data pnts_1000
read { z 3 }
limits 0 1 0 1
box

do i = 0,20 {
  data extract.$i
  read { c 1 }
  connect z c
}


