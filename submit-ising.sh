#!/bin/bash

set -ev

cargo build --release

STANDARD_FLAGS=( --movie-time '10^(1/4)' --translation-scale 1 --save-time 0.5 --max-iter 1e12 )

rq run -R -J ising-15-sad target/release/ising-mc --N 15 --sad-min-T 0.1 --save-as ising-15-sad.yaml ${STANDARD_FLAGS[@]}

rq run -R -J ising-15-wl target/release/ising-mc --N 15 --wl --save-as ising-15-wl.yaml ${STANDARD_FLAGS[@]}

rq run -R -J ising-15-samc-1e5 target/release/ising-mc --N 15 --samc-t0 1e5 --save-as ising-15-samc-1e5.yaml ${STANDARD_FLAGS[@]}

rq run -R -J ising-15-samc-1e4 target/release/ising-mc --N 15 --samc-t0 1e4 --save-as ising-15-samc-1e4.yaml ${STANDARD_FLAGS[@]}

rq run -R -J ising-15-samc-1e6 target/release/ising-mc --N 15 --samc-t0 1e6 --save-as ising-15-samc-1e6.yaml ${STANDARD_FLAGS[@]}
