#!/bin/bash

sbatch <<EOT
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1

RESULTS_DIRECTORY=results-mult$1-Nev$2-norm$3-nL$4

mkdir -p \${RESULTS_DIRECTORY}

./shifter.e RNG_mult=$1 RNG_nLoops=$2 shifter_norm=$3 shifter_nLoops=$4 \
            2> \${RESULTS_DIRECTORY}/run.out | grep ^OUT | awk '{print $4, $5}' \
            1> \${RESULTS_DIRECTORY}/pairs.out

exit 0
EOT
