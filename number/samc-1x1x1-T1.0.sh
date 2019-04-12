#!/bin/sh
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output samc-1x1x1-T1.0-t0-1e6.out

set -ev

hostname
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e3 --save-as samc-1x1x1-T1.0-t0-1e3.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5 > samc-1x1x1-T1.0-t0-1e3.out &
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e4 --save-as samc-1x1x1-T1.0-t0-1e4.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5 > samc-1x1x1-T1.0-t0-1e4.out &
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e5 --save-as samc-1x1x1-T1.0-t0-1e5.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5 > samc-1x1x1-T1.0-t0-1e5.out &


nice -19 ../target/release/number-mc --T 1 --samc-t0 1e6 --save-as samc-1x1x1-T1.0-t0-1e6.yaml --acceptance-rate 0.5 --cell-volume '4*1^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5 --save-time 0.5
