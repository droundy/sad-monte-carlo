#!/bin/sh
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output samc-v27-T1.0-t0-1e3.out

set -ev

hostname
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e3 --save-as samc-v27-T1.0-t0-1e3.yaml --acceptance-rate 0.5 --cell-volume 27 --movie-time "10^(1/4)" --well-width 1.3
