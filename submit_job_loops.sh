#!/bin/bash
#-------------------------------------------------------------------------------

#mult$1-Nev$2-norm$3-nL$4

for nL in 10 100 300 1000
do
  for mult in 3 4 5 7 10 25 50 75 100
  do
    for Nev in 100000 250000 500000
    do
      for norm in 0.1 0.25 0.5 0.75 1.0
      do
        bash submit_job.sh $mult $Nev $norm $nL
      done
    done
  done
done
