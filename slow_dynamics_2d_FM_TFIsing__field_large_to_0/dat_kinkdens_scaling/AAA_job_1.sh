#!/bin/bash

python=python
#prog=slow_dynamics_2d_FM_TFIsing.py
prog=../slow_dynamics_2d_FM_TFIsing.py

#----

date

for L in \
3 4
do
  for tau in \
   1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0 10.0 \
  11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 \
  30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0
  do
    echo "L=${L} tau=${tau}"
    ${python} ${prog} -Lx ${L} -Ly ${L} -tau ${tau} > dat_L${L}_tau${tau}
    mkdir -p fig_L${L}_tau${tau}
    mv fig_*.png fig_L${L}_tau${tau}
    date
  done
done
