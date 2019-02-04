#!/bin/sh
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output samc-t0-1e3-1x1x1.out

set -ev

hostname
nice -19 ../target/release/energy-number-mc --samc-t0 1e3 --save-as samc-t0-1e3-1x1x1.yaml --movie-name samc-t0-1e3-1x1x1 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --addremove-probability 0.2 --translation-scale 0.05
