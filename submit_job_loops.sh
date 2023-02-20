#!/bin/bash
#-------------------------------------------------------------------------------

#mult$1-Nev$2-norm$3-nL$4

for nL in 10 30 100 300 1000
do
  for mult in 3 4 5 7 10 15 20 25 30 40 50 60 75 100
  do
    for Nev in 100000 250000 500000 750000 1000000
    do
      for norm in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
      do
        bash submit_job.sh $mult $Nev $norm $nL
      done
    done
  done
done
