#!/bin/bash

for N in \
4 6 8 10 12 14 16
do
  grep '^kink density' dat_N${N}_tau* | sed 's/.*KinkDens//g' | sort -g -k 2 > dat_kinkdens_N${N}
done

gnuplot plot_png_kinkdens_vs_inversetau
gnuplot plot_png_kinkdens_vs_tau_loglog
