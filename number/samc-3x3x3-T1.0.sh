#!/bin/sh
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output samc-3x3x3-T1.0-t0-1e12.out

set -ev

hostname
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e10 --save-as samc-3x3x3-T1.0-t0-1e10.yaml --acceptance-rate 0.5 --cell-volume '4*3^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 > samc-3x3x3-T1.0-t0-1e10.out &
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e9 --save-as samc-3x3x3-T1.0-t0-1e9.yaml --acceptance-rate 0.5 --cell-volume '4*3^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 > samc-3x3x3-T1.0-t0-1e9.out &
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e11 --save-as samc-3x3x3-T1.0-t0-1e11.yaml --acceptance-rate 0.5 --cell-volume '4*3^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 > samc-3x3x3-T1.0-t0-1e11.out &

nice -19 ../target/release/number-mc --T 1 --samc-t0 1e12 --save-as samc-3x3x3-T1.0-t0-1e12.yaml --acceptance-rate 0.5 --cell-volume '4*3^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5
