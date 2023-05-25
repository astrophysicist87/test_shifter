#!/bin/bash

SHIFT_MODE=$1
MULTIPLICITY=$2
NUMBER_OF_EVENTS=$3
NUMBER_OF_LOOPS=$4

RESULTS_DIRECTORY=results-shiftmode_$SHIFT_MODE-mult$MULTIPLICITY-Nev$NUMBER_OF_EVENTS-nL$NUMBER_OF_LOOPS
rm -rf $RESULTS_DIRECTORY && mkdir -p $RESULTS_DIRECTORY

export OMP_NUM_THREADS=1

echo "Executing command: ./shifter.e $RESULTS_DIRECTORY $SHIFT_MODE RNG_mult=$MULTIPLICITY \
            RNG_xDir=1 RNG_yDir=1 RNG_nLoops=$NUMBER_OF_EVENTS \
            shifter_nLoops=$NUMBER_OF_LOOPS hybrid_cutoff=1000.0"

#valgrind \
./shifter.e $RESULTS_DIRECTORY $SHIFT_MODE RNG_mult=$MULTIPLICITY \
            RNG_xDir=1 RNG_yDir=1 RNG_nLoops=$NUMBER_OF_EVENTS \
            shifter_nLoops=$NUMBER_OF_LOOPS hybrid_cutoff=1000.0 #\
            #> $RESULTS_DIRECTORY/run.out
