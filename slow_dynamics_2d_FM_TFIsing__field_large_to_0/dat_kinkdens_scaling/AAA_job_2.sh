#!/bin/bash

for L in \
3 4
do
  grep '^kink density' dat_L${L}_tau* | sed 's/.*KinkDens//g' | sort -g -k 2 > dat_kinkdens_L${L}
done
