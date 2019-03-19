#!/bin/sh
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output samc-1x1x1-T1.0-t0-1e7.out

set -ev

hostname
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e10 --save-as samc-1x1x1-T1.0-t0-1e10.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5 > samc-1x1x1-T1.0-t0-1e10.out &
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e9 --save-as samc-1x1x1-T1.0-t0-1e9.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5 > samc-1x1x1-T1.0-t0-1e9.out &
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e8 --save-as samc-1x1x1-T1.0-t0-1e8.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5 > samc-1x1x1-T1.0-t0-1e8.out &


nice -19 ../target/release/number-mc --T 1 --samc-t0 1e7 --save-as samc-1x1x1-T1.0-t0-1e7.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5
