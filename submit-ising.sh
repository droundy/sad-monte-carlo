#!/bin/sh

set -ev

cargo build --release

STANDARD_FLAGS=" --movie-time '10^(1/4)' --translation-scale 1 --save-time 0.5"

rq run -R -J ising_15_sad target/release/ising-mc --N 15 --sad-min-T 0.1 --save-as ising_15_sad.yaml $STANDARD_FLAGS

rq run -R -J ising_15_wl target/release/ising-mc --N 15 --wl --save-as ising_15_wl.yaml $STANDARD_FLAGS
