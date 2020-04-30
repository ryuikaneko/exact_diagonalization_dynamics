#!/bin/bash

python=python
#prog=slow_dynamics_1d_FM_TFIsing.py
prog=../slow_dynamics_1d_FM_TFIsing.py

#----

date

for N in \
4 6 8 10 12 14 16
do
  for tau in \
   1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0 10.0 \
  11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 \
  30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0
  do
    echo "N=${N} tau=${tau}"
    ${python} ${prog} -N ${N} -tau ${tau} > dat_N${N}_tau${tau}
    mkdir -p fig_N${N}_tau${tau}
    mv fig_*.png fig_N${N}_tau${tau}
    date
  done
done
