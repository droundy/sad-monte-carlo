#!/bin/sh
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output samc-V40-T1.0-t0-1e9.out

set -ev

hostname
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e9 --save-as samc-V40-T1.0-t0-1e9.yaml --acceptance-rate 0.5 --cell-volume '4*(2.1834)^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5
nice -19 ../target/release/number-mc --T 1 --samc-t0 1e8 --save-as samc-V40-T1.0-t0-1e8.yaml --acceptance-rate 0.5 --cell-volume '4*(2.1834)^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5
nice -19 ../target/release/number-mc --T .9 --samc-t0 1e9 --save-as samc-V40-T0.9-t0-1e9.yaml --acceptance-rate 0.5 --cell-volume '4*(2.1834)^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5
nice -19 ../target/release/number-mc --T .9 --samc-t0 1e8 --save-as samc-V40-T0.9-t0-1e8.yaml --acceptance-rate 0.5 --cell-volume '4*(2.1834)^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5
nice -19 ../target/release/number-mc --T 1.1 --samc-t0 1e9 --save-as samc-V40-T1.1-t0-1e9.yaml --acceptance-rate 0.5 --cell-volume '4*(2.1834)^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5
nice -19 ../target/release/number-mc --T 1.1 --samc-t0 1e8 --save-as samc-V40-T1.1-t0-1e8.yaml --acceptance-rate 0.5 --cell-volume '4*(2.1834)^3*pi/6/0.545' --movie-time "10^(1/4)" --well-width 1.5
