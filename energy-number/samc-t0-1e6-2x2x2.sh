#!/bin/sh
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output samc-t0-1e6-2x2x2.out

set -ev

hostname
nice -19 ../target/release/energy-number-mc --samc-t0 1e6 --save-as samc-t0-1e6-2x2x2.yaml --addremove-probability 0.2  --cell-volume 30.74341 --movie-time "10^(1/4)" --well-width 1.3 --movie-name samc-t0-1e6-2x2x2 --save-time 0.5 --translation-scale 0.05
